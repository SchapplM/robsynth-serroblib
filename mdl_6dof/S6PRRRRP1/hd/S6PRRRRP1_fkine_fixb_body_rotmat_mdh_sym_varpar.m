% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRRRRP1 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:17
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6PRRRRP1_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP1_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRRP1_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP1_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:17:55
	% EndTime: 2020-11-04 21:17:55
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:17:55
	% EndTime: 2020-11-04 21:17:55
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t83 = cos(pkin(11));
	t82 = sin(pkin(11));
	t1 = [t83, -t82, 0, 0; t82, t83, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:17:55
	% EndTime: 2020-11-04 21:17:55
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t84 = sin(pkin(11));
	t85 = sin(pkin(6));
	t93 = t84 * t85;
	t86 = cos(pkin(11));
	t92 = t86 * t85;
	t87 = cos(pkin(6));
	t88 = sin(qJ(2));
	t91 = t87 * t88;
	t89 = cos(qJ(2));
	t90 = t87 * t89;
	t1 = [-t84 * t91 + t86 * t89, -t84 * t90 - t86 * t88, t93, t86 * pkin(1) + pkin(7) * t93 + 0; t84 * t89 + t86 * t91, -t84 * t88 + t86 * t90, -t92, t84 * pkin(1) - pkin(7) * t92 + 0; t85 * t88, t85 * t89, t87, t87 * pkin(7) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:17:55
	% EndTime: 2020-11-04 21:17:56
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (29->27), mult. (64->46), div. (0->0), fcn. (85->8), ass. (0->20)
	t97 = sin(pkin(6));
	t112 = t97 * pkin(7);
	t96 = sin(pkin(11));
	t99 = cos(pkin(6));
	t111 = t96 * t99;
	t98 = cos(pkin(11));
	t110 = t98 * t99;
	t100 = sin(qJ(3));
	t109 = t100 * t97;
	t102 = cos(qJ(3));
	t108 = t102 * t97;
	t101 = sin(qJ(2));
	t107 = t96 * t101;
	t103 = cos(qJ(2));
	t106 = t96 * t103;
	t105 = t98 * t101;
	t104 = t98 * t103;
	t95 = t99 * t105 + t106;
	t94 = t99 * t107 - t104;
	t1 = [-t94 * t102 + t96 * t109, t94 * t100 + t96 * t108, t99 * t106 + t105, (t98 * pkin(2) + pkin(8) * t111) * t103 + (-pkin(2) * t111 + t98 * pkin(8)) * t101 + t96 * t112 + t98 * pkin(1) + 0; t95 * t102 - t98 * t109, -t95 * t100 - t98 * t108, -t99 * t104 + t107, (t96 * pkin(2) - pkin(8) * t110) * t103 + (pkin(2) * t110 + t96 * pkin(8)) * t101 - t98 * t112 + t96 * pkin(1) + 0; t99 * t100 + t101 * t108, -t101 * t109 + t99 * t102, -t97 * t103, t99 * pkin(7) + qJ(1) + 0 + (pkin(2) * t101 - pkin(8) * t103) * t97; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:17:56
	% EndTime: 2020-11-04 21:17:56
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (50->27), mult. (75->40), div. (0->0), fcn. (100->10), ass. (0->22)
	t122 = sin(pkin(11));
	t123 = sin(pkin(6));
	t135 = t122 * t123;
	t124 = cos(pkin(11));
	t134 = t123 * t124;
	t127 = sin(qJ(2));
	t133 = t123 * t127;
	t125 = cos(pkin(6));
	t132 = t125 * t127;
	t128 = cos(qJ(2));
	t131 = t125 * t128;
	t130 = pkin(3) * sin(qJ(3)) + pkin(7);
	t129 = pkin(9) + pkin(8);
	t121 = qJ(3) + qJ(4);
	t120 = cos(t121);
	t119 = sin(t121);
	t118 = cos(qJ(3)) * pkin(3) + pkin(2);
	t116 = t122 * t131 + t124 * t127;
	t115 = t122 * t128 + t124 * t132;
	t114 = t122 * t127 - t124 * t131;
	t113 = t122 * t132 - t124 * t128;
	t1 = [-t113 * t120 + t119 * t135, t113 * t119 + t120 * t135, t116, t124 * pkin(1) - t113 * t118 + t116 * t129 + t130 * t135 + 0; t115 * t120 - t119 * t134, -t115 * t119 - t120 * t134, t114, t122 * pkin(1) + t114 * t129 + t115 * t118 - t130 * t134 + 0; t125 * t119 + t120 * t133, -t119 * t133 + t125 * t120, -t123 * t128, qJ(1) + 0 + t130 * t125 + (t118 * t127 - t128 * t129) * t123; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:17:56
	% EndTime: 2020-11-04 21:17:56
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (96->39), mult. (145->58), div. (0->0), fcn. (193->12), ass. (0->31)
	t151 = sin(pkin(11));
	t152 = sin(pkin(6));
	t167 = t151 * t152;
	t153 = cos(pkin(11));
	t166 = t152 * t153;
	t157 = sin(qJ(2));
	t165 = t152 * t157;
	t159 = cos(qJ(2));
	t164 = t152 * t159;
	t154 = cos(pkin(6));
	t163 = t154 * t157;
	t162 = t154 * t159;
	t161 = pkin(3) * sin(qJ(3)) + pkin(7);
	t160 = pkin(9) + pkin(8);
	t158 = cos(qJ(5));
	t155 = sin(qJ(5));
	t150 = qJ(3) + qJ(4);
	t149 = cos(t150);
	t148 = sin(t150);
	t147 = cos(qJ(3)) * pkin(3) + pkin(2);
	t145 = t151 * t162 + t153 * t157;
	t144 = t151 * t159 + t153 * t163;
	t143 = t151 * t157 - t153 * t162;
	t142 = t151 * t163 - t153 * t159;
	t141 = t154 * t148 + t149 * t165;
	t140 = t148 * t165 - t154 * t149;
	t139 = -t142 * t149 + t148 * t167;
	t138 = t144 * t149 - t148 * t166;
	t137 = t144 * t148 + t149 * t166;
	t136 = t142 * t148 + t149 * t167;
	t1 = [t139 * t158 + t145 * t155, -t139 * t155 + t145 * t158, -t136, t153 * pkin(1) + t139 * pkin(4) - t136 * pkin(10) - t142 * t147 + t145 * t160 + t161 * t167 + 0; t138 * t158 + t143 * t155, -t138 * t155 + t143 * t158, t137, t151 * pkin(1) + t138 * pkin(4) + t137 * pkin(10) + t143 * t160 + t144 * t147 - t161 * t166 + 0; t141 * t158 - t155 * t164, -t141 * t155 - t158 * t164, t140, t141 * pkin(4) + t140 * pkin(10) + qJ(1) + 0 + t161 * t154 + (t147 * t157 - t159 * t160) * t152; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:17:56
	% EndTime: 2020-11-04 21:17:56
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (107->42), mult. (161->60), div. (0->0), fcn. (211->12), ass. (0->33)
	t184 = sin(pkin(11));
	t185 = sin(pkin(6));
	t202 = t184 * t185;
	t186 = cos(pkin(11));
	t201 = t185 * t186;
	t191 = sin(qJ(2));
	t200 = t185 * t191;
	t193 = cos(qJ(2));
	t199 = t185 * t193;
	t187 = cos(pkin(6));
	t198 = t187 * t191;
	t197 = t187 * t193;
	t196 = pkin(3) * sin(qJ(3)) + pkin(7);
	t189 = sin(qJ(5));
	t195 = pkin(5) * t189 + pkin(8) + pkin(9);
	t192 = cos(qJ(5));
	t188 = -qJ(6) - pkin(10);
	t183 = qJ(3) + qJ(4);
	t182 = cos(t183);
	t181 = sin(t183);
	t180 = cos(qJ(3)) * pkin(3) + pkin(2);
	t179 = t192 * pkin(5) + pkin(4);
	t177 = t184 * t197 + t186 * t191;
	t176 = t184 * t193 + t186 * t198;
	t175 = t184 * t191 - t186 * t197;
	t174 = t184 * t198 - t186 * t193;
	t173 = t187 * t181 + t182 * t200;
	t172 = t181 * t200 - t187 * t182;
	t171 = -t174 * t182 + t181 * t202;
	t170 = t176 * t182 - t181 * t201;
	t169 = t176 * t181 + t182 * t201;
	t168 = t174 * t181 + t182 * t202;
	t1 = [t171 * t192 + t177 * t189, -t171 * t189 + t177 * t192, -t168, t186 * pkin(1) + t168 * t188 + t171 * t179 - t174 * t180 + t195 * t177 + t196 * t202 + 0; t170 * t192 + t175 * t189, -t170 * t189 + t175 * t192, t169, t184 * pkin(1) - t169 * t188 + t170 * t179 + t195 * t175 + t176 * t180 - t196 * t201 + 0; t173 * t192 - t189 * t199, -t173 * t189 - t192 * t199, t172, -t172 * t188 + t173 * t179 + qJ(1) + 0 + t196 * t187 + (t180 * t191 - t195 * t193) * t185; 0, 0, 0, 1;];
	Tc_mdh = t1;
end