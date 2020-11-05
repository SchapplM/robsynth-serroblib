% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRPRRR6 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:05
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6PRPRRR6_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR6_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRRR6_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR6_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:05:35
	% EndTime: 2020-11-04 21:05:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:05:35
	% EndTime: 2020-11-04 21:05:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t100 = cos(pkin(11));
	t99 = sin(pkin(11));
	t1 = [t100, -t99, 0, 0; t99, t100, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:05:35
	% EndTime: 2020-11-04 21:05:35
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t101 = sin(pkin(11));
	t102 = sin(pkin(6));
	t110 = t101 * t102;
	t103 = cos(pkin(11));
	t109 = t103 * t102;
	t104 = cos(pkin(6));
	t105 = sin(qJ(2));
	t108 = t104 * t105;
	t106 = cos(qJ(2));
	t107 = t104 * t106;
	t1 = [-t101 * t108 + t103 * t106, -t101 * t107 - t103 * t105, t110, t103 * pkin(1) + pkin(7) * t110 + 0; t101 * t106 + t103 * t108, -t101 * t105 + t103 * t107, -t109, t101 * pkin(1) - pkin(7) * t109 + 0; t102 * t105, t102 * t106, t104, t104 * pkin(7) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:05:35
	% EndTime: 2020-11-04 21:05:35
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (23->23), mult. (43->36), div. (0->0), fcn. (56->6), ass. (0->13)
	t111 = sin(pkin(11));
	t112 = sin(pkin(6));
	t122 = t111 * t112;
	t114 = cos(pkin(6));
	t121 = t111 * t114;
	t113 = cos(pkin(11));
	t120 = t113 * t112;
	t119 = t113 * t114;
	t115 = sin(qJ(2));
	t118 = t114 * t115;
	t116 = cos(qJ(2));
	t117 = t114 * t116;
	t1 = [t122, t111 * t118 - t113 * t116, t111 * t117 + t113 * t115, (t113 * pkin(2) + qJ(3) * t121) * t116 + (-pkin(2) * t121 + t113 * qJ(3)) * t115 + pkin(7) * t122 + t113 * pkin(1) + 0; -t120, -t111 * t116 - t113 * t118, t111 * t115 - t113 * t117, (t111 * pkin(2) - qJ(3) * t119) * t116 + (pkin(2) * t119 + t111 * qJ(3)) * t115 - pkin(7) * t120 + t111 * pkin(1) + 0; t114, -t112 * t115, -t112 * t116, t114 * pkin(7) + qJ(1) + 0 + (pkin(2) * t115 - qJ(3) * t116) * t112; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:05:35
	% EndTime: 2020-11-04 21:05:35
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (36->28), mult. (64->48), div. (0->0), fcn. (85->8), ass. (0->22)
	t125 = sin(pkin(11));
	t143 = t125 * qJ(3);
	t126 = sin(pkin(6));
	t129 = sin(qJ(4));
	t142 = t126 * t129;
	t131 = cos(qJ(4));
	t141 = t126 * t131;
	t132 = cos(qJ(2));
	t140 = t126 * t132;
	t133 = pkin(3) + pkin(7);
	t139 = t126 * t133;
	t127 = cos(pkin(11));
	t138 = t127 * qJ(3);
	t128 = cos(pkin(6));
	t130 = sin(qJ(2));
	t137 = t128 * t130;
	t136 = t128 * t132;
	t134 = pkin(2) + pkin(8);
	t135 = t128 * t134;
	t124 = t125 * t136 + t127 * t130;
	t123 = t125 * t130 - t127 * t136;
	t1 = [t124 * t129 + t125 * t141, t124 * t131 - t125 * t142, -t125 * t137 + t127 * t132, (t127 * t134 + t128 * t143) * t132 + (-t125 * t135 + t138) * t130 + t125 * t139 + t127 * pkin(1) + 0; t123 * t129 - t127 * t141, t123 * t131 + t127 * t142, t125 * t132 + t127 * t137, (t125 * t134 - t128 * t138) * t132 + (t127 * t135 + t143) * t130 - t127 * t139 + t125 * pkin(1) + 0; t128 * t131 - t129 * t140, -t128 * t129 - t131 * t140, t126 * t130, t133 * t128 + qJ(1) + 0 + (-qJ(3) * t132 + t130 * t134) * t126; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:05:35
	% EndTime: 2020-11-04 21:05:35
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (65->34), mult. (140->56), div. (0->0), fcn. (174->10), ass. (0->30)
	t150 = sin(pkin(6));
	t152 = cos(pkin(6));
	t155 = sin(qJ(2));
	t158 = cos(qJ(2));
	t160 = pkin(2) + pkin(8);
	t154 = sin(qJ(4));
	t157 = cos(qJ(4));
	t165 = t154 * pkin(4) - pkin(9) * t157 + qJ(3);
	t162 = -t160 * t155 + t165 * t158;
	t164 = t157 * pkin(4) + t154 * pkin(9) + pkin(3) + pkin(7);
	t176 = t164 * t150 + t162 * t152;
	t173 = t150 * t154;
	t172 = t150 * t155;
	t171 = t152 * t155;
	t170 = t152 * t158;
	t169 = t154 * t155;
	t168 = t154 * t158;
	t167 = t157 * t150;
	t163 = t152 * t168 + t167;
	t161 = t165 * t155 + t160 * t158 + pkin(1);
	t156 = cos(qJ(5));
	t153 = sin(qJ(5));
	t151 = cos(pkin(11));
	t149 = sin(pkin(11));
	t148 = -t150 * t168 + t152 * t157;
	t147 = t149 * t171 - t151 * t158;
	t146 = t149 * t158 + t151 * t171;
	t145 = -t149 * t169 + t163 * t151;
	t144 = t163 * t149 + t151 * t169;
	t1 = [t144 * t156 - t153 * t147, -t144 * t153 - t156 * t147, t149 * t173 - (t149 * t170 + t151 * t155) * t157, t176 * t149 + t161 * t151 + 0; -t145 * t156 + t146 * t153, t145 * t153 + t146 * t156, -t151 * t173 + (-t149 * t155 + t151 * t170) * t157, t161 * t149 - t176 * t151 + 0; t148 * t156 + t153 * t172, -t148 * t153 + t156 * t172, t152 * t154 + t158 * t167, -t162 * t150 + t164 * t152 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:05:35
	% EndTime: 2020-11-04 21:05:35
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (88->47), mult. (150->70), div. (0->0), fcn. (196->12), ass. (0->37)
	t213 = pkin(5) * sin(qJ(5));
	t191 = sin(pkin(11));
	t212 = t191 * qJ(3);
	t192 = sin(pkin(6));
	t196 = sin(qJ(4));
	t211 = t192 * t196;
	t197 = sin(qJ(2));
	t210 = t192 * t197;
	t198 = cos(qJ(4));
	t209 = t192 * t198;
	t199 = cos(qJ(2));
	t208 = t192 * t199;
	t201 = pkin(3) + pkin(7);
	t207 = t192 * t201;
	t193 = cos(pkin(11));
	t206 = t193 * qJ(3);
	t194 = cos(pkin(6));
	t205 = t194 * t197;
	t204 = t194 * t199;
	t202 = pkin(2) + pkin(8);
	t203 = t194 * t202;
	t200 = -pkin(10) - pkin(9);
	t190 = qJ(5) + qJ(6);
	t189 = cos(t190);
	t188 = sin(t190);
	t187 = cos(qJ(5)) * pkin(5) + pkin(4);
	t186 = t194 * t198 - t196 * t208;
	t185 = t194 * t196 + t198 * t208;
	t184 = -t191 * t205 + t193 * t199;
	t183 = t191 * t204 + t193 * t197;
	t182 = t191 * t199 + t193 * t205;
	t181 = t191 * t197 - t193 * t204;
	t180 = t181 * t196 - t193 * t209;
	t179 = t181 * t198 + t193 * t211;
	t178 = t183 * t196 + t191 * t209;
	t177 = -t183 * t198 + t191 * t211;
	t1 = [t178 * t189 + t184 * t188, -t178 * t188 + t184 * t189, t177, t178 * t187 - t177 * t200 + t184 * t213 + (t193 * t202 + t194 * t212) * t199 + (-t191 * t203 + t206) * t197 + t191 * t207 + t193 * pkin(1) + 0; t180 * t189 + t182 * t188, -t180 * t188 + t182 * t189, -t179, t180 * t187 + t179 * t200 + t182 * t213 + (t191 * t202 - t194 * t206) * t199 + (t193 * t203 + t212) * t197 - t193 * t207 + t191 * pkin(1) + 0; t186 * t189 + t188 * t210, -t186 * t188 + t189 * t210, t185, -t185 * t200 + t186 * t187 + t201 * t194 + qJ(1) + 0 + (-qJ(3) * t199 + (t202 + t213) * t197) * t192; 0, 0, 0, 1;];
	Tc_mdh = t1;
end