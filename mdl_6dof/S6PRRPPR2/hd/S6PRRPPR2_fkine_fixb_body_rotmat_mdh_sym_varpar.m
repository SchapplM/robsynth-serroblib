% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRRPPR2 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:06
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6PRRPPR2_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR2_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPPR2_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR2_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:06:52
	% EndTime: 2020-11-04 21:06:52
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:06:52
	% EndTime: 2020-11-04 21:06:52
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t89 = cos(pkin(10));
	t88 = sin(pkin(10));
	t1 = [t89, -t88, 0, 0; t88, t89, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:06:52
	% EndTime: 2020-11-04 21:06:52
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t90 = sin(pkin(10));
	t91 = sin(pkin(6));
	t99 = t90 * t91;
	t92 = cos(pkin(10));
	t98 = t92 * t91;
	t93 = cos(pkin(6));
	t94 = sin(qJ(2));
	t97 = t93 * t94;
	t95 = cos(qJ(2));
	t96 = t93 * t95;
	t1 = [-t90 * t97 + t92 * t95, -t90 * t96 - t92 * t94, t99, t92 * pkin(1) + pkin(7) * t99 + 0; t90 * t95 + t92 * t97, -t90 * t94 + t92 * t96, -t98, t90 * pkin(1) - pkin(7) * t98 + 0; t91 * t94, t91 * t95, t93, t93 * pkin(7) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:06:52
	% EndTime: 2020-11-04 21:06:53
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (29->27), mult. (64->48), div. (0->0), fcn. (85->8), ass. (0->18)
	t103 = sin(pkin(6));
	t116 = t103 * pkin(7);
	t102 = sin(pkin(10));
	t105 = cos(pkin(6));
	t115 = t102 * t105;
	t106 = sin(qJ(3));
	t114 = t103 * t106;
	t108 = cos(qJ(3));
	t113 = t103 * t108;
	t104 = cos(pkin(10));
	t112 = t104 * t105;
	t107 = sin(qJ(2));
	t111 = t105 * t107;
	t109 = cos(qJ(2));
	t110 = t105 * t109;
	t101 = t102 * t109 + t104 * t111;
	t100 = t102 * t111 - t104 * t109;
	t1 = [-t100 * t108 + t102 * t114, t100 * t106 + t102 * t113, t102 * t110 + t104 * t107, (t104 * pkin(2) + pkin(8) * t115) * t109 + (-pkin(2) * t115 + t104 * pkin(8)) * t107 + t102 * t116 + t104 * pkin(1) + 0; t101 * t108 - t104 * t114, -t101 * t106 - t104 * t113, t102 * t107 - t104 * t110, (t102 * pkin(2) - pkin(8) * t112) * t109 + (pkin(2) * t112 + t102 * pkin(8)) * t107 - t104 * t116 + t102 * pkin(1) + 0; t105 * t106 + t107 * t113, t105 * t108 - t107 * t114, -t103 * t109, t105 * pkin(7) + qJ(1) + 0 + (pkin(2) * t107 - pkin(8) * t109) * t103; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:06:53
	% EndTime: 2020-11-04 21:06:53
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (50->27), mult. (75->40), div. (0->0), fcn. (100->10), ass. (0->22)
	t126 = sin(pkin(10));
	t127 = sin(pkin(6));
	t139 = t126 * t127;
	t128 = cos(pkin(10));
	t138 = t127 * t128;
	t132 = sin(qJ(2));
	t137 = t127 * t132;
	t129 = cos(pkin(6));
	t136 = t129 * t132;
	t133 = cos(qJ(2));
	t135 = t129 * t133;
	t134 = pkin(3) * sin(qJ(3)) + pkin(7);
	t130 = qJ(4) + pkin(8);
	t125 = qJ(3) + pkin(11);
	t124 = cos(t125);
	t123 = sin(t125);
	t122 = cos(qJ(3)) * pkin(3) + pkin(2);
	t120 = t126 * t135 + t128 * t132;
	t119 = t126 * t133 + t128 * t136;
	t118 = t126 * t132 - t128 * t135;
	t117 = t126 * t136 - t128 * t133;
	t1 = [-t117 * t124 + t123 * t139, t117 * t123 + t124 * t139, t120, t128 * pkin(1) - t117 * t122 + t120 * t130 + t134 * t139 + 0; t119 * t124 - t123 * t138, -t119 * t123 - t124 * t138, t118, t126 * pkin(1) + t118 * t130 + t119 * t122 - t134 * t138 + 0; t129 * t123 + t124 * t137, -t123 * t137 + t129 * t124, -t127 * t133, qJ(1) + 0 + t134 * t129 + (t122 * t132 - t130 * t133) * t127; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:06:53
	% EndTime: 2020-11-04 21:06:53
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (78->37), mult. (111->46), div. (0->0), fcn. (146->10), ass. (0->28)
	t155 = sin(pkin(10));
	t156 = sin(pkin(6));
	t168 = t155 * t156;
	t157 = cos(pkin(10));
	t167 = t156 * t157;
	t161 = sin(qJ(2));
	t166 = t156 * t161;
	t158 = cos(pkin(6));
	t165 = t158 * t161;
	t162 = cos(qJ(2));
	t164 = t158 * t162;
	t163 = pkin(3) * sin(qJ(3)) + pkin(7);
	t159 = qJ(4) + pkin(8);
	t154 = qJ(3) + pkin(11);
	t153 = cos(t154);
	t152 = sin(t154);
	t151 = cos(qJ(3)) * pkin(3) + pkin(2);
	t149 = t155 * t164 + t157 * t161;
	t148 = t155 * t162 + t157 * t165;
	t147 = t155 * t161 - t157 * t164;
	t146 = t155 * t165 - t157 * t162;
	t145 = t158 * t152 + t153 * t166;
	t144 = t152 * t166 - t158 * t153;
	t143 = -t146 * t153 + t152 * t168;
	t142 = t148 * t153 - t152 * t167;
	t141 = t148 * t152 + t153 * t167;
	t140 = t146 * t152 + t153 * t168;
	t1 = [t149, -t143, -t140, t157 * pkin(1) + t143 * pkin(4) - t140 * qJ(5) - t146 * t151 + t149 * t159 + t163 * t168 + 0; t147, -t142, t141, t155 * pkin(1) + t142 * pkin(4) + t141 * qJ(5) + t147 * t159 + t148 * t151 - t163 * t167 + 0; -t156 * t162, -t145, t144, t145 * pkin(4) + t144 * qJ(5) + qJ(1) + 0 + t163 * t158 + (t151 * t161 - t159 * t162) * t156; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:06:53
	% EndTime: 2020-11-04 21:06:53
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (125->40), mult. (170->62), div. (0->0), fcn. (204->14), ass. (0->39)
	t181 = sin(pkin(11));
	t184 = cos(pkin(11));
	t193 = pkin(4) + pkin(9);
	t175 = qJ(5) * t181 + t193 * t184 + pkin(3);
	t176 = qJ(5) * t184 - t181 * t193;
	t183 = sin(pkin(6));
	t186 = cos(pkin(6));
	t189 = sin(qJ(2));
	t191 = cos(qJ(3));
	t188 = sin(qJ(3));
	t198 = t175 * t188 + pkin(7);
	t199 = t176 * t188 + pkin(2);
	t201 = t186 * t189;
	t179 = qJ(4) + pkin(5) + pkin(8);
	t192 = cos(qJ(2));
	t206 = t179 * t192;
	t210 = (t199 * t189 - t206) * t186 - t198 * t183 + (t175 * t201 + t183 * t176) * t191;
	t182 = sin(pkin(10));
	t205 = t182 * t183;
	t185 = cos(pkin(10));
	t204 = t183 * t185;
	t203 = t183 * t189;
	t202 = t183 * t192;
	t200 = t186 * t192;
	t197 = t175 * t191 + t199;
	t171 = t182 * t201 - t185 * t192;
	t180 = qJ(3) + pkin(11);
	t177 = sin(t180);
	t178 = cos(t180);
	t196 = -t171 * t177 - t178 * t205;
	t173 = t182 * t192 + t185 * t201;
	t195 = t173 * t177 + t178 * t204;
	t194 = t179 * t189 + t197 * t192 + pkin(1);
	t190 = cos(qJ(6));
	t187 = sin(qJ(6));
	t174 = t182 * t200 + t185 * t189;
	t172 = t182 * t189 - t185 * t200;
	t170 = t177 * t203 - t186 * t178;
	t1 = [t174 * t190 + t196 * t187, -t174 * t187 + t196 * t190, -t171 * t178 + t177 * t205, -t210 * t182 + t194 * t185 + 0; t172 * t190 + t195 * t187, -t172 * t187 + t195 * t190, t173 * t178 - t177 * t204, t194 * t182 + t210 * t185 + 0; t170 * t187 - t190 * t202, t170 * t190 + t187 * t202, t186 * t177 + t178 * t203, qJ(1) + 0 + (t197 * t189 - t206) * t183 + (-t176 * t191 + t198) * t186; 0, 0, 0, 1;];
	Tc_mdh = t1;
end