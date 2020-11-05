% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRRPRR2 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:10
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6PRRPRR2_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR2_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRR2_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR2_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [12x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:10:40
	% EndTime: 2020-11-04 21:10:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:10:41
	% EndTime: 2020-11-04 21:10:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t82 = cos(pkin(11));
	t81 = sin(pkin(11));
	t1 = [t82, -t81, 0, 0; t81, t82, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:10:41
	% EndTime: 2020-11-04 21:10:41
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t83 = sin(pkin(11));
	t84 = sin(pkin(6));
	t92 = t83 * t84;
	t85 = cos(pkin(11));
	t91 = t85 * t84;
	t86 = cos(pkin(6));
	t87 = sin(qJ(2));
	t90 = t86 * t87;
	t88 = cos(qJ(2));
	t89 = t86 * t88;
	t1 = [-t83 * t90 + t85 * t88, -t83 * t89 - t85 * t87, t92, t85 * pkin(1) + pkin(7) * t92 + 0; t83 * t88 + t85 * t90, -t83 * t87 + t85 * t89, -t91, t83 * pkin(1) - pkin(7) * t91 + 0; t84 * t87, t84 * t88, t86, t86 * pkin(7) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:10:41
	% EndTime: 2020-11-04 21:10:41
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (29->27), mult. (64->46), div. (0->0), fcn. (85->8), ass. (0->20)
	t96 = sin(pkin(6));
	t111 = t96 * pkin(7);
	t95 = sin(pkin(11));
	t98 = cos(pkin(6));
	t110 = t95 * t98;
	t99 = sin(qJ(3));
	t109 = t96 * t99;
	t97 = cos(pkin(11));
	t108 = t97 * t98;
	t101 = cos(qJ(3));
	t107 = t101 * t96;
	t100 = sin(qJ(2));
	t106 = t95 * t100;
	t102 = cos(qJ(2));
	t105 = t95 * t102;
	t104 = t97 * t100;
	t103 = t97 * t102;
	t94 = t98 * t104 + t105;
	t93 = t98 * t106 - t103;
	t1 = [-t93 * t101 + t95 * t109, t95 * t107 + t93 * t99, t98 * t105 + t104, (t97 * pkin(2) + pkin(8) * t110) * t102 + (-pkin(2) * t110 + t97 * pkin(8)) * t100 + t95 * t111 + t97 * pkin(1) + 0; t94 * t101 - t97 * t109, -t97 * t107 - t94 * t99, -t98 * t103 + t106, (t95 * pkin(2) - pkin(8) * t108) * t102 + (pkin(2) * t108 + t95 * pkin(8)) * t100 - t97 * t111 + t95 * pkin(1) + 0; t100 * t107 + t98 * t99, -t100 * t109 + t98 * t101, -t96 * t102, t98 * pkin(7) + qJ(1) + 0 + (pkin(2) * t100 - pkin(8) * t102) * t96; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:10:41
	% EndTime: 2020-11-04 21:10:41
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (50->27), mult. (75->40), div. (0->0), fcn. (100->10), ass. (0->22)
	t121 = sin(pkin(11));
	t122 = sin(pkin(6));
	t134 = t121 * t122;
	t123 = cos(pkin(11));
	t133 = t122 * t123;
	t127 = sin(qJ(2));
	t132 = t122 * t127;
	t124 = cos(pkin(6));
	t131 = t124 * t127;
	t128 = cos(qJ(2));
	t130 = t124 * t128;
	t129 = pkin(3) * sin(qJ(3)) + pkin(7);
	t125 = qJ(4) + pkin(8);
	t120 = qJ(3) + pkin(12);
	t119 = cos(t120);
	t118 = sin(t120);
	t117 = cos(qJ(3)) * pkin(3) + pkin(2);
	t115 = t121 * t130 + t123 * t127;
	t114 = t121 * t128 + t123 * t131;
	t113 = t121 * t127 - t123 * t130;
	t112 = t121 * t131 - t123 * t128;
	t1 = [-t112 * t119 + t118 * t134, t112 * t118 + t119 * t134, t115, t123 * pkin(1) - t112 * t117 + t115 * t125 + t129 * t134 + 0; t114 * t119 - t118 * t133, -t114 * t118 - t119 * t133, t113, t121 * pkin(1) + t113 * t125 + t114 * t117 - t129 * t133 + 0; t124 * t118 + t119 * t132, -t118 * t132 + t124 * t119, -t122 * t128, qJ(1) + 0 + t129 * t124 + (t117 * t127 - t125 * t128) * t122; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:10:41
	% EndTime: 2020-11-04 21:10:41
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (96->39), mult. (145->58), div. (0->0), fcn. (193->12), ass. (0->31)
	t150 = sin(pkin(11));
	t151 = sin(pkin(6));
	t166 = t150 * t151;
	t152 = cos(pkin(11));
	t165 = t151 * t152;
	t157 = sin(qJ(2));
	t164 = t151 * t157;
	t159 = cos(qJ(2));
	t163 = t151 * t159;
	t153 = cos(pkin(6));
	t162 = t153 * t157;
	t161 = t153 * t159;
	t160 = pkin(3) * sin(qJ(3)) + pkin(7);
	t158 = cos(qJ(5));
	t155 = sin(qJ(5));
	t154 = qJ(4) + pkin(8);
	t149 = qJ(3) + pkin(12);
	t148 = cos(t149);
	t147 = sin(t149);
	t146 = cos(qJ(3)) * pkin(3) + pkin(2);
	t144 = t150 * t161 + t152 * t157;
	t143 = t150 * t159 + t152 * t162;
	t142 = t150 * t157 - t152 * t161;
	t141 = t150 * t162 - t152 * t159;
	t140 = t153 * t147 + t148 * t164;
	t139 = t147 * t164 - t153 * t148;
	t138 = -t141 * t148 + t147 * t166;
	t137 = t143 * t148 - t147 * t165;
	t136 = t143 * t147 + t148 * t165;
	t135 = t141 * t147 + t148 * t166;
	t1 = [t138 * t158 + t144 * t155, -t138 * t155 + t144 * t158, -t135, t152 * pkin(1) + t138 * pkin(4) - t135 * pkin(9) - t141 * t146 + t144 * t154 + t160 * t166 + 0; t137 * t158 + t142 * t155, -t137 * t155 + t142 * t158, t136, t150 * pkin(1) + t137 * pkin(4) + t136 * pkin(9) + t142 * t154 + t143 * t146 - t160 * t165 + 0; t140 * t158 - t155 * t163, -t140 * t155 - t158 * t163, t139, t140 * pkin(4) + t139 * pkin(9) + qJ(1) + 0 + t160 * t153 + (t146 * t157 - t154 * t159) * t151; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:10:41
	% EndTime: 2020-11-04 21:10:41
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (119->43), mult. (161->60), div. (0->0), fcn. (211->14), ass. (0->34)
	t186 = sin(pkin(11));
	t187 = sin(pkin(6));
	t203 = t186 * t187;
	t188 = cos(pkin(11));
	t202 = t187 * t188;
	t193 = sin(qJ(2));
	t201 = t187 * t193;
	t194 = cos(qJ(2));
	t200 = t187 * t194;
	t189 = cos(pkin(6));
	t199 = t189 * t193;
	t198 = t189 * t194;
	t197 = pkin(3) * sin(qJ(3)) + pkin(7);
	t196 = pkin(5) * sin(qJ(5)) + pkin(8) + qJ(4);
	t195 = -pkin(10) - pkin(9);
	t185 = qJ(5) + qJ(6);
	t184 = qJ(3) + pkin(12);
	t183 = cos(t185);
	t182 = sin(t185);
	t181 = cos(t184);
	t180 = sin(t184);
	t179 = cos(qJ(3)) * pkin(3) + pkin(2);
	t178 = cos(qJ(5)) * pkin(5) + pkin(4);
	t176 = t186 * t198 + t188 * t193;
	t175 = t186 * t194 + t188 * t199;
	t174 = t186 * t193 - t188 * t198;
	t173 = t186 * t199 - t188 * t194;
	t172 = t189 * t180 + t181 * t201;
	t171 = t180 * t201 - t189 * t181;
	t170 = -t173 * t181 + t180 * t203;
	t169 = t175 * t181 - t180 * t202;
	t168 = t175 * t180 + t181 * t202;
	t167 = t173 * t180 + t181 * t203;
	t1 = [t170 * t183 + t176 * t182, -t170 * t182 + t176 * t183, -t167, t188 * pkin(1) + t167 * t195 + t170 * t178 - t173 * t179 + t196 * t176 + t197 * t203 + 0; t169 * t183 + t174 * t182, -t169 * t182 + t174 * t183, t168, t186 * pkin(1) - t168 * t195 + t169 * t178 + t196 * t174 + t175 * t179 - t197 * t202 + 0; t172 * t183 - t182 * t200, -t172 * t182 - t183 * t200, t171, -t171 * t195 + t172 * t178 + qJ(1) + 0 + t197 * t189 + (t179 * t193 - t196 * t194) * t187; 0, 0, 0, 1;];
	Tc_mdh = t1;
end