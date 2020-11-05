% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRRRR5 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:47
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRRRRR5_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR5_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRR5_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR5_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [12x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:47:20
	% EndTime: 2020-11-04 22:47:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:47:20
	% EndTime: 2020-11-04 22:47:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t91 = cos(qJ(1));
	t90 = sin(qJ(1));
	t1 = [t91, -t90, 0, 0; t90, t91, 0, 0; 0, 0, 1, pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:47:20
	% EndTime: 2020-11-04 22:47:20
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t92 = sin(pkin(6));
	t95 = sin(qJ(1));
	t103 = t95 * t92;
	t94 = sin(qJ(2));
	t102 = t95 * t94;
	t96 = cos(qJ(2));
	t101 = t95 * t96;
	t97 = cos(qJ(1));
	t100 = t97 * t92;
	t99 = t97 * t94;
	t98 = t97 * t96;
	t93 = cos(pkin(6));
	t1 = [-t93 * t102 + t98, -t93 * t101 - t99, t103, t97 * pkin(1) + pkin(8) * t103 + 0; t93 * t99 + t101, t93 * t98 - t102, -t100, t95 * pkin(1) - pkin(8) * t100 + 0; t92 * t94, t92 * t96, t93, t93 * pkin(8) + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:47:20
	% EndTime: 2020-11-04 22:47:20
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (29->23), mult. (60->37), div. (0->0), fcn. (81->8), ass. (0->20)
	t107 = sin(pkin(6));
	t112 = cos(qJ(3));
	t122 = t107 * t112;
	t109 = sin(qJ(3));
	t121 = t109 * t107;
	t110 = sin(qJ(2));
	t120 = t110 * t112;
	t111 = sin(qJ(1));
	t119 = t111 * t110;
	t113 = cos(qJ(2));
	t118 = t111 * t113;
	t114 = cos(qJ(1));
	t117 = t114 * t110;
	t116 = t114 * t113;
	t115 = pkin(2) * t110 - pkin(9) * t113;
	t108 = cos(pkin(6));
	t106 = t113 * pkin(2) + t110 * pkin(9) + pkin(1);
	t105 = -t108 * t117 - t118;
	t104 = t107 * pkin(8) - t115 * t108;
	t1 = [(-t108 * t120 + t121) * t111 + t112 * t116, (t108 * t119 - t116) * t109 + t111 * t122, t108 * t118 + t117, t104 * t111 + t106 * t114 + 0; -t105 * t112 - t114 * t121, t105 * t109 - t114 * t122, -t108 * t116 + t119, -t104 * t114 + t106 * t111 + 0; t107 * t120 + t108 * t109, t108 * t112 - t110 * t121, -t107 * t113, t108 * pkin(8) + t115 * t107 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:47:20
	% EndTime: 2020-11-04 22:47:21
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (54->26), mult. (69->38), div. (0->0), fcn. (90->10), ass. (0->25)
	t131 = sin(pkin(6));
	t134 = sin(qJ(2));
	t147 = t131 * t134;
	t135 = sin(qJ(1));
	t146 = t131 * t135;
	t137 = cos(qJ(1));
	t145 = t131 * t137;
	t144 = t135 * t134;
	t136 = cos(qJ(2));
	t143 = t135 * t136;
	t142 = t137 * t134;
	t141 = t137 * t136;
	t140 = sin(qJ(3)) * pkin(3) + pkin(8);
	t127 = cos(qJ(3)) * pkin(3) + pkin(2);
	t138 = pkin(10) + pkin(9);
	t139 = t127 * t134 - t138 * t136;
	t132 = cos(pkin(6));
	t130 = qJ(3) + qJ(4);
	t129 = cos(t130);
	t128 = sin(t130);
	t126 = t132 * t142 + t143;
	t125 = t132 * t144 - t141;
	t124 = t127 * t136 + t138 * t134 + pkin(1);
	t123 = t131 * t140 - t139 * t132;
	t1 = [-t125 * t129 + t128 * t146, t125 * t128 + t129 * t146, t132 * t143 + t142, t123 * t135 + t124 * t137 + 0; t126 * t129 - t128 * t145, -t126 * t128 - t129 * t145, -t132 * t141 + t144, -t123 * t137 + t124 * t135 + 0; t132 * t128 + t129 * t147, -t128 * t147 + t132 * t129, -t131 * t136, t139 * t131 + t140 * t132 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:47:21
	% EndTime: 2020-11-04 22:47:21
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (77->30), mult. (81->40), div. (0->0), fcn. (106->12), ass. (0->25)
	t158 = qJ(3) + qJ(4);
	t172 = pkin(8) + pkin(4) * sin(t158) + sin(qJ(3)) * pkin(3);
	t159 = sin(pkin(6));
	t161 = sin(qJ(2));
	t171 = t159 * t161;
	t162 = sin(qJ(1));
	t170 = t159 * t162;
	t164 = cos(qJ(1));
	t169 = t159 * t164;
	t168 = t162 * t161;
	t163 = cos(qJ(2));
	t167 = t162 * t163;
	t166 = t164 * t161;
	t165 = t164 * t163;
	t160 = cos(pkin(6));
	t157 = -pkin(11) - pkin(10) - pkin(9);
	t156 = qJ(5) + t158;
	t155 = cos(t156);
	t154 = sin(t156);
	t152 = pkin(4) * cos(t158) + cos(qJ(3)) * pkin(3) + pkin(2);
	t151 = -t160 * t168 + t165;
	t150 = t160 * t167 + t166;
	t149 = t160 * t166 + t167;
	t148 = -t160 * t165 + t168;
	t1 = [t151 * t155 + t154 * t170, -t151 * t154 + t155 * t170, t150, t164 * pkin(1) - t150 * t157 + t151 * t152 + t172 * t170 + 0; t149 * t155 - t154 * t169, -t149 * t154 - t155 * t169, t148, t162 * pkin(1) - t148 * t157 + t149 * t152 - t172 * t169 + 0; t160 * t154 + t155 * t171, -t154 * t171 + t160 * t155, -t159 * t163, pkin(7) + 0 + t172 * t160 + (t152 * t161 + t157 * t163) * t159; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:47:21
	% EndTime: 2020-11-04 22:47:21
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (141->41), mult. (151->58), div. (0->0), fcn. (199->14), ass. (0->34)
	t189 = qJ(3) + qJ(4);
	t206 = pkin(8) + pkin(4) * sin(t189) + sin(qJ(3)) * pkin(3);
	t190 = sin(pkin(6));
	t193 = sin(qJ(2));
	t205 = t190 * t193;
	t194 = sin(qJ(1));
	t204 = t190 * t194;
	t196 = cos(qJ(2));
	t203 = t190 * t196;
	t197 = cos(qJ(1));
	t202 = t190 * t197;
	t201 = t194 * t193;
	t200 = t194 * t196;
	t199 = t197 * t193;
	t198 = t197 * t196;
	t195 = cos(qJ(6));
	t192 = sin(qJ(6));
	t191 = cos(pkin(6));
	t188 = -pkin(11) - pkin(10) - pkin(9);
	t187 = qJ(5) + t189;
	t186 = cos(t187);
	t185 = sin(t187);
	t183 = pkin(4) * cos(t189) + cos(qJ(3)) * pkin(3) + pkin(2);
	t182 = -t191 * t201 + t198;
	t181 = t191 * t200 + t199;
	t180 = t191 * t199 + t200;
	t179 = -t191 * t198 + t201;
	t178 = t191 * t185 + t186 * t205;
	t177 = t185 * t205 - t191 * t186;
	t176 = t182 * t186 + t185 * t204;
	t175 = t182 * t185 - t186 * t204;
	t174 = t180 * t186 - t185 * t202;
	t173 = t180 * t185 + t186 * t202;
	t1 = [t176 * t195 + t181 * t192, -t176 * t192 + t181 * t195, t175, t197 * pkin(1) + t176 * pkin(5) + t175 * pkin(12) - t181 * t188 + t182 * t183 + t206 * t204 + 0; t174 * t195 + t179 * t192, -t174 * t192 + t179 * t195, t173, t194 * pkin(1) + t174 * pkin(5) + t173 * pkin(12) - t179 * t188 + t180 * t183 - t206 * t202 + 0; t178 * t195 - t192 * t203, -t178 * t192 - t195 * t203, t177, t178 * pkin(5) + t177 * pkin(12) + pkin(7) + 0 + t206 * t191 + (t183 * t193 + t188 * t196) * t190; 0, 0, 0, 1;];
	Tc_mdh = t1;
end