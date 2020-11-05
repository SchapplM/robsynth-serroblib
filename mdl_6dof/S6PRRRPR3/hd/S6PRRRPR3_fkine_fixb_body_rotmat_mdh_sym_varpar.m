% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRRRPR3 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:15
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6PRRRPR3_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR3_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRPR3_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR3_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:15:26
	% EndTime: 2020-11-04 21:15:26
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:15:26
	% EndTime: 2020-11-04 21:15:26
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t90 = cos(pkin(11));
	t89 = sin(pkin(11));
	t1 = [t90, -t89, 0, 0; t89, t90, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:15:26
	% EndTime: 2020-11-04 21:15:26
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t91 = sin(pkin(11));
	t92 = sin(pkin(6));
	t100 = t91 * t92;
	t93 = cos(pkin(11));
	t99 = t93 * t92;
	t94 = cos(pkin(6));
	t95 = sin(qJ(2));
	t98 = t94 * t95;
	t96 = cos(qJ(2));
	t97 = t94 * t96;
	t1 = [-t91 * t98 + t93 * t96, -t91 * t97 - t93 * t95, t100, t93 * pkin(1) + pkin(7) * t100 + 0; t91 * t96 + t93 * t98, -t91 * t95 + t93 * t97, -t99, t91 * pkin(1) - pkin(7) * t99 + 0; t92 * t95, t92 * t96, t94, t94 * pkin(7) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:15:26
	% EndTime: 2020-11-04 21:15:26
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (29->27), mult. (64->48), div. (0->0), fcn. (85->8), ass. (0->18)
	t104 = sin(pkin(6));
	t117 = t104 * pkin(7);
	t103 = sin(pkin(11));
	t106 = cos(pkin(6));
	t116 = t103 * t106;
	t107 = sin(qJ(3));
	t115 = t104 * t107;
	t109 = cos(qJ(3));
	t114 = t104 * t109;
	t105 = cos(pkin(11));
	t113 = t105 * t106;
	t108 = sin(qJ(2));
	t112 = t106 * t108;
	t110 = cos(qJ(2));
	t111 = t106 * t110;
	t102 = t103 * t110 + t105 * t112;
	t101 = t103 * t112 - t105 * t110;
	t1 = [-t101 * t109 + t103 * t115, t101 * t107 + t103 * t114, t103 * t111 + t105 * t108, (t105 * pkin(2) + pkin(8) * t116) * t110 + (-pkin(2) * t116 + t105 * pkin(8)) * t108 + t103 * t117 + t105 * pkin(1) + 0; t102 * t109 - t105 * t115, -t102 * t107 - t105 * t114, t103 * t108 - t105 * t111, (t103 * pkin(2) - pkin(8) * t113) * t110 + (pkin(2) * t113 + t103 * pkin(8)) * t108 - t105 * t117 + t103 * pkin(1) + 0; t106 * t107 + t108 * t114, t106 * t109 - t108 * t115, -t104 * t110, t106 * pkin(7) + qJ(1) + 0 + (pkin(2) * t108 - pkin(8) * t110) * t104; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:15:26
	% EndTime: 2020-11-04 21:15:26
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (50->27), mult. (75->40), div. (0->0), fcn. (100->10), ass. (0->22)
	t127 = sin(pkin(11));
	t128 = sin(pkin(6));
	t140 = t127 * t128;
	t129 = cos(pkin(11));
	t139 = t128 * t129;
	t132 = sin(qJ(2));
	t138 = t128 * t132;
	t130 = cos(pkin(6));
	t137 = t130 * t132;
	t133 = cos(qJ(2));
	t136 = t130 * t133;
	t135 = pkin(3) * sin(qJ(3)) + pkin(7);
	t134 = pkin(9) + pkin(8);
	t126 = qJ(3) + qJ(4);
	t125 = cos(t126);
	t124 = sin(t126);
	t123 = cos(qJ(3)) * pkin(3) + pkin(2);
	t121 = t127 * t136 + t129 * t132;
	t120 = t127 * t133 + t129 * t137;
	t119 = t127 * t132 - t129 * t136;
	t118 = t127 * t137 - t129 * t133;
	t1 = [-t118 * t125 + t124 * t140, t118 * t124 + t125 * t140, t121, t129 * pkin(1) - t118 * t123 + t121 * t134 + t135 * t140 + 0; t120 * t125 - t124 * t139, -t120 * t124 - t125 * t139, t119, t127 * pkin(1) + t119 * t134 + t120 * t123 - t135 * t139 + 0; t130 * t124 + t125 * t138, -t124 * t138 + t130 * t125, -t128 * t133, qJ(1) + 0 + t135 * t130 + (t123 * t132 - t133 * t134) * t128; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:15:26
	% EndTime: 2020-11-04 21:15:26
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (78->37), mult. (111->46), div. (0->0), fcn. (146->10), ass. (0->28)
	t156 = sin(pkin(11));
	t157 = sin(pkin(6));
	t169 = t156 * t157;
	t158 = cos(pkin(11));
	t168 = t157 * t158;
	t161 = sin(qJ(2));
	t167 = t157 * t161;
	t159 = cos(pkin(6));
	t166 = t159 * t161;
	t162 = cos(qJ(2));
	t165 = t159 * t162;
	t164 = pkin(3) * sin(qJ(3)) + pkin(7);
	t163 = pkin(9) + pkin(8);
	t155 = qJ(3) + qJ(4);
	t154 = cos(t155);
	t153 = sin(t155);
	t152 = cos(qJ(3)) * pkin(3) + pkin(2);
	t150 = t156 * t165 + t158 * t161;
	t149 = t156 * t162 + t158 * t166;
	t148 = t156 * t161 - t158 * t165;
	t147 = t156 * t166 - t158 * t162;
	t146 = t159 * t153 + t154 * t167;
	t145 = t153 * t167 - t159 * t154;
	t144 = -t147 * t154 + t153 * t169;
	t143 = t149 * t154 - t153 * t168;
	t142 = t149 * t153 + t154 * t168;
	t141 = t147 * t153 + t154 * t169;
	t1 = [t150, -t144, -t141, t158 * pkin(1) + t144 * pkin(4) - t141 * qJ(5) - t147 * t152 + t150 * t163 + t164 * t169 + 0; t148, -t143, t142, t156 * pkin(1) + t143 * pkin(4) + t142 * qJ(5) + t148 * t163 + t149 * t152 - t164 * t168 + 0; -t157 * t162, -t146, t145, t146 * pkin(4) + t145 * qJ(5) + qJ(1) + 0 + t164 * t159 + (t152 * t161 - t162 * t163) * t157; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:15:26
	% EndTime: 2020-11-04 21:15:26
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (125->40), mult. (170->62), div. (0->0), fcn. (204->14), ass. (0->39)
	t191 = cos(qJ(4));
	t208 = sin(qJ(4));
	t209 = pkin(4) + pkin(10);
	t176 = qJ(5) * t208 + t209 * t191 + pkin(3);
	t177 = qJ(5) * t191 - t208 * t209;
	t184 = sin(pkin(6));
	t186 = cos(pkin(6));
	t189 = sin(qJ(2));
	t192 = cos(qJ(3));
	t188 = sin(qJ(3));
	t198 = t176 * t188 + pkin(7);
	t199 = t177 * t188 + pkin(2);
	t201 = t186 * t189;
	t181 = pkin(5) + pkin(8) + pkin(9);
	t193 = cos(qJ(2));
	t206 = t181 * t193;
	t212 = (t199 * t189 - t206) * t186 - t198 * t184 + (t176 * t201 + t177 * t184) * t192;
	t183 = sin(pkin(11));
	t205 = t183 * t184;
	t185 = cos(pkin(11));
	t204 = t184 * t185;
	t203 = t184 * t189;
	t202 = t184 * t193;
	t200 = t186 * t193;
	t197 = t176 * t192 + t199;
	t172 = t183 * t201 - t185 * t193;
	t182 = qJ(3) + qJ(4);
	t179 = sin(t182);
	t180 = cos(t182);
	t196 = -t172 * t179 - t180 * t205;
	t174 = t183 * t193 + t185 * t201;
	t195 = t174 * t179 + t180 * t204;
	t194 = t181 * t189 + t197 * t193 + pkin(1);
	t190 = cos(qJ(6));
	t187 = sin(qJ(6));
	t175 = t183 * t200 + t185 * t189;
	t173 = t183 * t189 - t185 * t200;
	t171 = t179 * t203 - t186 * t180;
	t1 = [t175 * t190 + t196 * t187, -t175 * t187 + t196 * t190, -t172 * t180 + t179 * t205, -t212 * t183 + t194 * t185 + 0; t173 * t190 + t195 * t187, -t173 * t187 + t195 * t190, t174 * t180 - t179 * t204, t194 * t183 + t212 * t185 + 0; t171 * t187 - t190 * t202, t171 * t190 + t187 * t202, t186 * t179 + t180 * t203, qJ(1) + 0 + (t197 * t189 - t206) * t184 + (-t177 * t192 + t198) * t186; 0, 0, 0, 1;];
	Tc_mdh = t1;
end