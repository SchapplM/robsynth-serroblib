% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRRRRP2 (for one body)
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
% Datum: 2020-11-04 21:18
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6PRRRRP2_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP2_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRRP2_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP2_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:18:17
	% EndTime: 2020-11-04 21:18:17
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:18:17
	% EndTime: 2020-11-04 21:18:17
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t85 = cos(pkin(11));
	t84 = sin(pkin(11));
	t1 = [t85, -t84, 0, 0; t84, t85, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:18:17
	% EndTime: 2020-11-04 21:18:17
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t86 = sin(pkin(11));
	t87 = sin(pkin(6));
	t95 = t86 * t87;
	t88 = cos(pkin(11));
	t94 = t88 * t87;
	t89 = cos(pkin(6));
	t90 = sin(qJ(2));
	t93 = t89 * t90;
	t91 = cos(qJ(2));
	t92 = t89 * t91;
	t1 = [-t86 * t93 + t88 * t91, -t86 * t92 - t88 * t90, t95, t88 * pkin(1) + pkin(7) * t95 + 0; t86 * t91 + t88 * t93, -t86 * t90 + t88 * t92, -t94, t86 * pkin(1) - pkin(7) * t94 + 0; t87 * t90, t87 * t91, t89, t89 * pkin(7) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:18:17
	% EndTime: 2020-11-04 21:18:17
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (29->27), mult. (64->48), div. (0->0), fcn. (85->8), ass. (0->18)
	t99 = sin(pkin(6));
	t112 = t99 * pkin(7);
	t101 = cos(pkin(6));
	t98 = sin(pkin(11));
	t111 = t101 * t98;
	t102 = sin(qJ(3));
	t110 = t102 * t99;
	t104 = cos(qJ(3));
	t109 = t104 * t99;
	t100 = cos(pkin(11));
	t108 = t100 * t101;
	t103 = sin(qJ(2));
	t107 = t101 * t103;
	t105 = cos(qJ(2));
	t106 = t101 * t105;
	t97 = t100 * t107 + t98 * t105;
	t96 = -t100 * t105 + t98 * t107;
	t1 = [-t96 * t104 + t98 * t110, t96 * t102 + t98 * t109, t100 * t103 + t98 * t106, (t100 * pkin(2) + pkin(8) * t111) * t105 + (-pkin(2) * t111 + t100 * pkin(8)) * t103 + t98 * t112 + t100 * pkin(1) + 0; -t100 * t110 + t97 * t104, -t100 * t109 - t97 * t102, -t100 * t106 + t98 * t103, (t98 * pkin(2) - pkin(8) * t108) * t105 + (pkin(2) * t108 + t98 * pkin(8)) * t103 - t100 * t112 + t98 * pkin(1) + 0; t101 * t102 + t103 * t109, t101 * t104 - t103 * t110, -t99 * t105, t101 * pkin(7) + qJ(1) + 0 + (pkin(2) * t103 - pkin(8) * t105) * t99; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:18:17
	% EndTime: 2020-11-04 21:18:17
	% DurationCPUTime: 0.12s
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
	% StartTime: 2020-11-04 21:18:17
	% EndTime: 2020-11-04 21:18:17
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
	% StartTime: 2020-11-04 21:18:17
	% EndTime: 2020-11-04 21:18:17
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (134->45), mult. (207->64), div. (0->0), fcn. (275->12), ass. (0->37)
	t189 = sin(pkin(11));
	t190 = sin(pkin(6));
	t205 = t189 * t190;
	t191 = cos(pkin(11));
	t204 = t190 * t191;
	t195 = sin(qJ(2));
	t203 = t190 * t195;
	t197 = cos(qJ(2));
	t202 = t190 * t197;
	t192 = cos(pkin(6));
	t201 = t192 * t195;
	t200 = t192 * t197;
	t199 = pkin(3) * sin(qJ(3)) + pkin(7);
	t198 = pkin(9) + pkin(8);
	t196 = cos(qJ(5));
	t193 = sin(qJ(5));
	t188 = qJ(3) + qJ(4);
	t187 = cos(t188);
	t186 = sin(t188);
	t185 = cos(qJ(3)) * pkin(3) + pkin(2);
	t183 = t189 * t200 + t191 * t195;
	t182 = t189 * t197 + t191 * t201;
	t181 = t189 * t195 - t191 * t200;
	t180 = t189 * t201 - t191 * t197;
	t179 = t192 * t186 + t187 * t203;
	t178 = t186 * t203 - t192 * t187;
	t177 = t179 * t196 - t193 * t202;
	t176 = t179 * t193 + t196 * t202;
	t175 = -t180 * t187 + t186 * t205;
	t174 = t182 * t187 - t186 * t204;
	t173 = t182 * t186 + t187 * t204;
	t172 = t180 * t186 + t187 * t205;
	t171 = t175 * t196 + t183 * t193;
	t170 = t175 * t193 - t183 * t196;
	t169 = t174 * t196 + t181 * t193;
	t168 = t174 * t193 - t181 * t196;
	t1 = [t171, -t172, t170, t191 * pkin(1) + t175 * pkin(4) + t171 * pkin(5) - t172 * pkin(10) + t170 * qJ(6) - t180 * t185 + t183 * t198 + t199 * t205 + 0; t169, t173, t168, t189 * pkin(1) + t174 * pkin(4) + t169 * pkin(5) + t173 * pkin(10) + t168 * qJ(6) + t181 * t198 + t182 * t185 - t199 * t204 + 0; t177, t178, t176, t179 * pkin(4) + t177 * pkin(5) + t178 * pkin(10) + t176 * qJ(6) + qJ(1) + 0 + t199 * t192 + (t185 * t195 - t197 * t198) * t190; 0, 0, 0, 1;];
	Tc_mdh = t1;
end