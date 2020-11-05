% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRPRRP5 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:02
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6PRPRRP5_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP5_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRRP5_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP5_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:02:56
	% EndTime: 2020-11-04 21:02:56
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:02:56
	% EndTime: 2020-11-04 21:02:56
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t99 = cos(pkin(10));
	t98 = sin(pkin(10));
	t1 = [t99, -t98, 0, 0; t98, t99, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:02:56
	% EndTime: 2020-11-04 21:02:56
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t100 = sin(pkin(10));
	t101 = sin(pkin(6));
	t109 = t100 * t101;
	t102 = cos(pkin(10));
	t108 = t102 * t101;
	t103 = cos(pkin(6));
	t104 = sin(qJ(2));
	t107 = t103 * t104;
	t105 = cos(qJ(2));
	t106 = t103 * t105;
	t1 = [-t100 * t107 + t102 * t105, -t100 * t106 - t102 * t104, t109, t102 * pkin(1) + pkin(7) * t109 + 0; t100 * t105 + t102 * t107, -t100 * t104 + t102 * t106, -t108, t100 * pkin(1) - pkin(7) * t108 + 0; t101 * t104, t101 * t105, t103, t103 * pkin(7) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:02:56
	% EndTime: 2020-11-04 21:02:56
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (23->23), mult. (43->36), div. (0->0), fcn. (56->6), ass. (0->13)
	t110 = sin(pkin(10));
	t111 = sin(pkin(6));
	t121 = t110 * t111;
	t113 = cos(pkin(6));
	t120 = t110 * t113;
	t112 = cos(pkin(10));
	t119 = t112 * t111;
	t118 = t112 * t113;
	t114 = sin(qJ(2));
	t117 = t113 * t114;
	t115 = cos(qJ(2));
	t116 = t113 * t115;
	t1 = [t121, t110 * t117 - t112 * t115, t110 * t116 + t112 * t114, (t112 * pkin(2) + qJ(3) * t120) * t115 + (-pkin(2) * t120 + t112 * qJ(3)) * t114 + pkin(7) * t121 + t112 * pkin(1) + 0; -t119, -t110 * t115 - t112 * t117, t110 * t114 - t112 * t116, (t110 * pkin(2) - qJ(3) * t118) * t115 + (pkin(2) * t118 + t110 * qJ(3)) * t114 - pkin(7) * t119 + t110 * pkin(1) + 0; t113, -t111 * t114, -t111 * t115, t113 * pkin(7) + qJ(1) + 0 + (pkin(2) * t114 - qJ(3) * t115) * t111; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:02:56
	% EndTime: 2020-11-04 21:02:56
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (36->28), mult. (64->48), div. (0->0), fcn. (85->8), ass. (0->22)
	t124 = sin(pkin(10));
	t142 = t124 * qJ(3);
	t125 = sin(pkin(6));
	t128 = sin(qJ(4));
	t141 = t125 * t128;
	t130 = cos(qJ(4));
	t140 = t125 * t130;
	t131 = cos(qJ(2));
	t139 = t125 * t131;
	t132 = pkin(3) + pkin(7);
	t138 = t125 * t132;
	t126 = cos(pkin(10));
	t137 = t126 * qJ(3);
	t127 = cos(pkin(6));
	t129 = sin(qJ(2));
	t136 = t127 * t129;
	t135 = t127 * t131;
	t133 = pkin(2) + pkin(8);
	t134 = t127 * t133;
	t123 = t124 * t135 + t126 * t129;
	t122 = t124 * t129 - t126 * t135;
	t1 = [t123 * t128 + t124 * t140, t123 * t130 - t124 * t141, -t124 * t136 + t126 * t131, (t126 * t133 + t127 * t142) * t131 + (-t124 * t134 + t137) * t129 + t124 * t138 + t126 * pkin(1) + 0; t122 * t128 - t126 * t140, t122 * t130 + t126 * t141, t124 * t131 + t126 * t136, (t124 * t133 - t127 * t137) * t131 + (t126 * t134 + t142) * t129 - t126 * t138 + t124 * pkin(1) + 0; t127 * t130 - t128 * t139, -t127 * t128 - t130 * t139, t125 * t129, t132 * t127 + qJ(1) + 0 + (-qJ(3) * t131 + t129 * t133) * t125; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:02:56
	% EndTime: 2020-11-04 21:02:56
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (65->34), mult. (140->56), div. (0->0), fcn. (174->10), ass. (0->30)
	t149 = sin(pkin(6));
	t151 = cos(pkin(6));
	t154 = sin(qJ(2));
	t157 = cos(qJ(2));
	t159 = pkin(2) + pkin(8);
	t153 = sin(qJ(4));
	t156 = cos(qJ(4));
	t164 = t153 * pkin(4) - pkin(9) * t156 + qJ(3);
	t161 = -t159 * t154 + t164 * t157;
	t163 = t156 * pkin(4) + t153 * pkin(9) + pkin(3) + pkin(7);
	t175 = t163 * t149 + t161 * t151;
	t172 = t149 * t153;
	t171 = t149 * t154;
	t170 = t151 * t154;
	t169 = t151 * t157;
	t168 = t153 * t154;
	t167 = t153 * t157;
	t166 = t156 * t149;
	t162 = t151 * t167 + t166;
	t160 = t164 * t154 + t159 * t157 + pkin(1);
	t155 = cos(qJ(5));
	t152 = sin(qJ(5));
	t150 = cos(pkin(10));
	t148 = sin(pkin(10));
	t147 = -t149 * t167 + t151 * t156;
	t146 = t148 * t170 - t150 * t157;
	t145 = t148 * t157 + t150 * t170;
	t144 = -t148 * t168 + t162 * t150;
	t143 = t162 * t148 + t150 * t168;
	t1 = [t143 * t155 - t152 * t146, -t143 * t152 - t155 * t146, t148 * t172 - (t148 * t169 + t150 * t154) * t156, t175 * t148 + t160 * t150 + 0; -t144 * t155 + t145 * t152, t144 * t152 + t145 * t155, -t150 * t172 + (-t148 * t154 + t150 * t169) * t156, t160 * t148 - t175 * t150 + 0; t147 * t155 + t152 * t171, -t147 * t152 + t155 * t171, t151 * t153 + t157 * t166, -t161 * t149 + t163 * t151 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:02:56
	% EndTime: 2020-11-04 21:02:56
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (76->46), mult. (150->70), div. (0->0), fcn. (196->10), ass. (0->37)
	t187 = sin(pkin(10));
	t189 = cos(pkin(10));
	t197 = cos(qJ(2));
	t190 = cos(pkin(6));
	t194 = sin(qJ(2));
	t202 = t190 * t194;
	t181 = t187 * t197 + t189 * t202;
	t192 = sin(qJ(5));
	t211 = t181 * t192;
	t183 = -t187 * t202 + t189 * t197;
	t210 = t183 * t192;
	t209 = t187 * qJ(3);
	t188 = sin(pkin(6));
	t193 = sin(qJ(4));
	t208 = t188 * t193;
	t207 = t188 * t194;
	t196 = cos(qJ(4));
	t206 = t188 * t196;
	t205 = t188 * t197;
	t198 = pkin(3) + pkin(7);
	t204 = t188 * t198;
	t203 = t189 * qJ(3);
	t201 = t190 * t197;
	t199 = pkin(2) + pkin(8);
	t200 = t190 * t199;
	t195 = cos(qJ(5));
	t191 = -qJ(6) - pkin(9);
	t186 = t195 * pkin(5) + pkin(4);
	t185 = t190 * t196 - t193 * t205;
	t184 = t190 * t193 + t196 * t205;
	t182 = t187 * t201 + t189 * t194;
	t180 = t187 * t194 - t189 * t201;
	t179 = t180 * t193 - t189 * t206;
	t178 = t180 * t196 + t189 * t208;
	t177 = t182 * t193 + t187 * t206;
	t176 = -t182 * t196 + t187 * t208;
	t1 = [t177 * t195 + t210, -t177 * t192 + t183 * t195, t176, t177 * t186 - t176 * t191 + pkin(5) * t210 + (t189 * t199 + t190 * t209) * t197 + (-t187 * t200 + t203) * t194 + t187 * t204 + t189 * pkin(1) + 0; t179 * t195 + t211, -t179 * t192 + t181 * t195, -t178, t179 * t186 + t178 * t191 + pkin(5) * t211 + (t187 * t199 - t190 * t203) * t197 + (t189 * t200 + t209) * t194 - t189 * t204 + t187 * pkin(1) + 0; t185 * t195 + t192 * t207, -t185 * t192 + t195 * t207, t184, -t184 * t191 + t185 * t186 + t198 * t190 + qJ(1) + 0 + (-qJ(3) * t197 + (pkin(5) * t192 + t199) * t194) * t188; 0, 0, 0, 1;];
	Tc_mdh = t1;
end