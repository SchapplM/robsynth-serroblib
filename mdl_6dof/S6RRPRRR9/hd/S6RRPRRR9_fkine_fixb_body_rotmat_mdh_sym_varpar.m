% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPRRR9 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:19
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRPRRR9_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR9_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRR9_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR9_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [12x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:19:40
	% EndTime: 2020-11-04 22:19:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:19:40
	% EndTime: 2020-11-04 22:19:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t87 = cos(qJ(1));
	t86 = sin(qJ(1));
	t1 = [t87, -t86, 0, 0; t86, t87, 0, 0; 0, 0, 1, pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:19:40
	% EndTime: 2020-11-04 22:19:40
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t88 = sin(pkin(6));
	t91 = sin(qJ(1));
	t99 = t91 * t88;
	t90 = sin(qJ(2));
	t98 = t91 * t90;
	t92 = cos(qJ(2));
	t97 = t91 * t92;
	t93 = cos(qJ(1));
	t96 = t93 * t88;
	t95 = t93 * t90;
	t94 = t93 * t92;
	t89 = cos(pkin(6));
	t1 = [-t89 * t98 + t94, -t89 * t97 - t95, t99, t93 * pkin(1) + pkin(8) * t99 + 0; t89 * t95 + t97, t89 * t94 - t98, -t96, t91 * pkin(1) - pkin(8) * t96 + 0; t88 * t90, t88 * t92, t89, t89 * pkin(8) + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:19:40
	% EndTime: 2020-11-04 22:19:40
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (29->22), mult. (60->36), div. (0->0), fcn. (81->8), ass. (0->21)
	t105 = sin(pkin(6));
	t108 = sin(qJ(2));
	t119 = t105 * t108;
	t109 = sin(qJ(1));
	t118 = t105 * t109;
	t111 = cos(qJ(1));
	t117 = t105 * t111;
	t116 = t109 * t108;
	t110 = cos(qJ(2));
	t115 = t109 * t110;
	t114 = t111 * t108;
	t113 = t111 * t110;
	t112 = pkin(2) * t108 - qJ(3) * t110;
	t107 = cos(pkin(6));
	t106 = cos(pkin(12));
	t104 = sin(pkin(12));
	t103 = t110 * pkin(2) + t108 * qJ(3) + pkin(1);
	t102 = t107 * t114 + t115;
	t101 = t107 * t116 - t113;
	t100 = t105 * pkin(8) - t112 * t107;
	t1 = [-t101 * t106 + t104 * t118, t101 * t104 + t106 * t118, t107 * t115 + t114, t100 * t109 + t103 * t111 + 0; t102 * t106 - t104 * t117, -t102 * t104 - t106 * t117, -t107 * t113 + t116, -t100 * t111 + t103 * t109 + 0; t107 * t104 + t106 * t119, -t104 * t119 + t107 * t106, -t105 * t110, t107 * pkin(8) + t112 * t105 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:19:40
	% EndTime: 2020-11-04 22:19:40
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (54->26), mult. (68->38), div. (0->0), fcn. (89->10), ass. (0->25)
	t129 = sin(pkin(6));
	t132 = sin(qJ(2));
	t143 = t129 * t132;
	t133 = sin(qJ(1));
	t142 = t129 * t133;
	t135 = cos(qJ(1));
	t141 = t129 * t135;
	t140 = t133 * t132;
	t134 = cos(qJ(2));
	t139 = t133 * t134;
	t138 = t135 * t132;
	t137 = t135 * t134;
	t125 = cos(pkin(12)) * pkin(3) + pkin(2);
	t131 = qJ(3) + pkin(9);
	t136 = t125 * t132 - t131 * t134;
	t130 = cos(pkin(6));
	t128 = pkin(12) + qJ(4);
	t127 = cos(t128);
	t126 = sin(t128);
	t124 = sin(pkin(12)) * pkin(3) + pkin(8);
	t123 = t130 * t138 + t139;
	t122 = t130 * t140 - t137;
	t121 = t125 * t134 + t131 * t132 + pkin(1);
	t120 = t129 * t124 - t136 * t130;
	t1 = [-t122 * t127 + t126 * t142, t122 * t126 + t127 * t142, t130 * t139 + t138, t120 * t133 + t121 * t135 + 0; t123 * t127 - t126 * t141, -t123 * t126 - t127 * t141, -t130 * t137 + t140, -t120 * t135 + t121 * t133 + 0; t130 * t126 + t127 * t143, -t126 * t143 + t130 * t127, -t129 * t134, t124 * t130 + t136 * t129 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:19:40
	% EndTime: 2020-11-04 22:19:40
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (77->30), mult. (81->40), div. (0->0), fcn. (106->12), ass. (0->25)
	t154 = pkin(12) + qJ(4);
	t168 = pkin(8) + pkin(4) * sin(t154) + sin(pkin(12)) * pkin(3);
	t155 = sin(pkin(6));
	t157 = sin(qJ(2));
	t167 = t155 * t157;
	t158 = sin(qJ(1));
	t166 = t155 * t158;
	t160 = cos(qJ(1));
	t165 = t155 * t160;
	t164 = t158 * t157;
	t159 = cos(qJ(2));
	t163 = t158 * t159;
	t162 = t160 * t157;
	t161 = t160 * t159;
	t156 = cos(pkin(6));
	t153 = -pkin(10) - pkin(9) - qJ(3);
	t152 = qJ(5) + t154;
	t151 = cos(t152);
	t150 = sin(t152);
	t148 = pkin(4) * cos(t154) + cos(pkin(12)) * pkin(3) + pkin(2);
	t147 = -t156 * t164 + t161;
	t146 = t156 * t163 + t162;
	t145 = t156 * t162 + t163;
	t144 = -t156 * t161 + t164;
	t1 = [t147 * t151 + t150 * t166, -t147 * t150 + t151 * t166, t146, t160 * pkin(1) - t146 * t153 + t147 * t148 + t168 * t166 + 0; t145 * t151 - t150 * t165, -t145 * t150 - t151 * t165, t144, t158 * pkin(1) - t144 * t153 + t145 * t148 - t168 * t165 + 0; t156 * t150 + t151 * t167, -t150 * t167 + t156 * t151, -t155 * t159, pkin(7) + 0 + t168 * t156 + (t148 * t157 + t153 * t159) * t155; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:19:40
	% EndTime: 2020-11-04 22:19:40
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (141->41), mult. (151->58), div. (0->0), fcn. (199->14), ass. (0->34)
	t185 = pkin(12) + qJ(4);
	t202 = pkin(8) + pkin(4) * sin(t185) + sin(pkin(12)) * pkin(3);
	t186 = sin(pkin(6));
	t189 = sin(qJ(2));
	t201 = t186 * t189;
	t190 = sin(qJ(1));
	t200 = t186 * t190;
	t192 = cos(qJ(2));
	t199 = t186 * t192;
	t193 = cos(qJ(1));
	t198 = t186 * t193;
	t197 = t190 * t189;
	t196 = t190 * t192;
	t195 = t193 * t189;
	t194 = t193 * t192;
	t191 = cos(qJ(6));
	t188 = sin(qJ(6));
	t187 = cos(pkin(6));
	t184 = -pkin(10) - pkin(9) - qJ(3);
	t183 = qJ(5) + t185;
	t182 = cos(t183);
	t181 = sin(t183);
	t179 = pkin(4) * cos(t185) + cos(pkin(12)) * pkin(3) + pkin(2);
	t178 = -t187 * t197 + t194;
	t177 = t187 * t196 + t195;
	t176 = t187 * t195 + t196;
	t175 = -t187 * t194 + t197;
	t174 = t187 * t181 + t182 * t201;
	t173 = t181 * t201 - t187 * t182;
	t172 = t178 * t182 + t181 * t200;
	t171 = t178 * t181 - t182 * t200;
	t170 = t176 * t182 - t181 * t198;
	t169 = t176 * t181 + t182 * t198;
	t1 = [t172 * t191 + t177 * t188, -t172 * t188 + t177 * t191, t171, t193 * pkin(1) + t172 * pkin(5) + t171 * pkin(11) - t177 * t184 + t178 * t179 + t202 * t200 + 0; t170 * t191 + t175 * t188, -t170 * t188 + t175 * t191, t169, t190 * pkin(1) + t170 * pkin(5) + t169 * pkin(11) - t175 * t184 + t176 * t179 - t202 * t198 + 0; t174 * t191 - t188 * t199, -t174 * t188 - t191 * t199, t173, t174 * pkin(5) + t173 * pkin(11) + pkin(7) + 0 + t202 * t187 + (t179 * t189 + t184 * t192) * t186; 0, 0, 0, 1;];
	Tc_mdh = t1;
end