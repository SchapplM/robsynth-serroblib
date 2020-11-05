% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRPPR10 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:25
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRRPPR10_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR10_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPPR10_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR10_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:25:38
	% EndTime: 2020-11-04 22:25:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:25:38
	% EndTime: 2020-11-04 22:25:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t95 = cos(qJ(1));
	t94 = sin(qJ(1));
	t1 = [t95, -t94, 0, 0; t94, t95, 0, 0; 0, 0, 1, pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:25:38
	% EndTime: 2020-11-04 22:25:38
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t96 = sin(pkin(6));
	t99 = sin(qJ(1));
	t107 = t99 * t96;
	t98 = sin(qJ(2));
	t106 = t99 * t98;
	t101 = cos(qJ(1));
	t105 = t101 * t96;
	t104 = t101 * t98;
	t100 = cos(qJ(2));
	t103 = t99 * t100;
	t102 = t101 * t100;
	t97 = cos(pkin(6));
	t1 = [-t97 * t106 + t102, -t97 * t103 - t104, t107, t101 * pkin(1) + pkin(8) * t107 + 0; t97 * t104 + t103, t97 * t102 - t106, -t105, t99 * pkin(1) - pkin(8) * t105 + 0; t96 * t98, t96 * t100, t97, t97 * pkin(8) + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:25:38
	% EndTime: 2020-11-04 22:25:38
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (29->23), mult. (60->37), div. (0->0), fcn. (81->8), ass. (0->20)
	t111 = sin(pkin(6));
	t116 = cos(qJ(3));
	t126 = t111 * t116;
	t113 = sin(qJ(3));
	t125 = t113 * t111;
	t114 = sin(qJ(2));
	t124 = t114 * t116;
	t115 = sin(qJ(1));
	t123 = t115 * t114;
	t117 = cos(qJ(2));
	t122 = t115 * t117;
	t118 = cos(qJ(1));
	t121 = t118 * t114;
	t120 = t118 * t117;
	t119 = pkin(2) * t114 - pkin(9) * t117;
	t112 = cos(pkin(6));
	t110 = t117 * pkin(2) + t114 * pkin(9) + pkin(1);
	t109 = -t112 * t121 - t122;
	t108 = t111 * pkin(8) - t119 * t112;
	t1 = [(-t112 * t124 + t125) * t115 + t116 * t120, (t112 * t123 - t120) * t113 + t115 * t126, t112 * t122 + t121, t108 * t115 + t110 * t118 + 0; -t109 * t116 - t118 * t125, t109 * t113 - t118 * t126, -t112 * t120 + t123, -t108 * t118 + t110 * t115 + 0; t111 * t124 + t112 * t113, t112 * t116 - t114 * t125, -t111 * t117, t112 * pkin(8) + t119 * t111 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:25:38
	% EndTime: 2020-11-04 22:25:38
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (45->27), mult. (78->41), div. (0->0), fcn. (99->8), ass. (0->22)
	t131 = sin(pkin(6));
	t136 = cos(qJ(3));
	t147 = t131 * t136;
	t133 = sin(qJ(3));
	t146 = t133 * t131;
	t134 = sin(qJ(2));
	t145 = t134 * t136;
	t135 = sin(qJ(1));
	t144 = t135 * t134;
	t137 = cos(qJ(2));
	t143 = t135 * t137;
	t138 = cos(qJ(1));
	t142 = t138 * t134;
	t141 = t138 * t137;
	t130 = t136 * pkin(3) + qJ(4) * t133 + pkin(2);
	t140 = t137 * pkin(9) - t130 * t134;
	t139 = t133 * pkin(3) - qJ(4) * t136 + pkin(8);
	t132 = cos(pkin(6));
	t129 = -t132 * t142 - t143;
	t128 = t134 * pkin(9) + t130 * t137 + pkin(1);
	t127 = t131 * t139 + t140 * t132;
	t1 = [t132 * t143 + t142, (t132 * t145 - t146) * t135 - t136 * t141, (-t132 * t144 + t141) * t133 - t135 * t147, t127 * t135 + t128 * t138 + 0; -t132 * t141 + t144, t129 * t136 + t138 * t146, -t129 * t133 + t138 * t147, -t127 * t138 + t128 * t135 + 0; -t131 * t137, -t131 * t145 - t132 * t133, -t132 * t136 + t134 * t146, -t140 * t131 + t139 * t132 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:25:38
	% EndTime: 2020-11-04 22:25:38
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (70->34), mult. (112->56), div. (0->0), fcn. (146->10), ass. (0->31)
	t158 = sin(pkin(6));
	t166 = cos(qJ(2));
	t177 = t158 * t166;
	t160 = cos(pkin(6));
	t165 = cos(qJ(3));
	t176 = t160 * t165;
	t175 = t160 * t166;
	t162 = sin(qJ(3));
	t174 = t162 * t158;
	t163 = sin(qJ(2));
	t173 = t162 * t163;
	t172 = t162 * t166;
	t171 = t165 * t158;
	t161 = qJ(5) + pkin(3);
	t156 = qJ(4) * t162 + t161 * t165 + pkin(2);
	t168 = pkin(4) + pkin(9);
	t170 = t156 * t163 - t168 * t166;
	t169 = qJ(4) * t165 - t161 * t162 - pkin(8);
	t167 = cos(qJ(1));
	t164 = sin(qJ(1));
	t159 = cos(pkin(11));
	t157 = sin(pkin(11));
	t155 = t158 * t173 - t176;
	t154 = t160 * t173 + t171;
	t153 = t157 * t163 - t159 * t172;
	t152 = t157 * t172 + t159 * t163;
	t151 = t156 * t166 + t168 * t163 + pkin(1);
	t150 = -t157 * t154 + t159 * t175;
	t149 = t159 * t154 + t157 * t175;
	t148 = -t169 * t158 - t170 * t160;
	t1 = [t150 * t164 + t167 * t152, -t149 * t164 - t167 * t153, (-t163 * t176 + t174) * t164 + t165 * t166 * t167, t148 * t164 + t151 * t167 + 0; -t150 * t167 + t164 * t152, t149 * t167 - t164 * t153, (t167 * t160 * t163 + t164 * t166) * t165 - t167 * t174, -t148 * t167 + t151 * t164 + 0; t155 * t157 - t159 * t177, t155 * t159 + t157 * t177, t160 * t162 + t163 * t171, t170 * t158 - t169 * t160 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:25:38
	% EndTime: 2020-11-04 22:25:38
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (103->37), mult. (125->55), div. (0->0), fcn. (159->12), ass. (0->34)
	t193 = sin(pkin(6));
	t199 = cos(qJ(2));
	t210 = t193 * t199;
	t194 = cos(pkin(6));
	t198 = cos(qJ(3));
	t209 = t194 * t198;
	t195 = sin(qJ(3));
	t208 = t195 * t193;
	t196 = sin(qJ(2));
	t207 = t195 * t196;
	t197 = sin(qJ(1));
	t206 = t197 * t199;
	t205 = t198 * t193;
	t200 = cos(qJ(1));
	t204 = t199 * t200;
	t203 = t200 * t196;
	t188 = sin(pkin(11)) * pkin(5) + qJ(4);
	t191 = qJ(5) + pkin(3) + pkin(10);
	t182 = t188 * t195 + t191 * t198 + pkin(2);
	t187 = cos(pkin(11)) * pkin(5) + pkin(4) + pkin(9);
	t202 = t182 * t196 - t187 * t199;
	t201 = t188 * t198 - t191 * t195 - pkin(8);
	t192 = pkin(11) + qJ(6);
	t190 = cos(t192);
	t189 = sin(t192);
	t186 = t194 * t206 + t203;
	t185 = t194 * t204 - t197 * t196;
	t184 = t193 * t207 - t209;
	t183 = t194 * t207 + t205;
	t181 = -t183 * t197 + t195 * t204;
	t180 = t183 * t200 + t195 * t206;
	t179 = t182 * t199 + t187 * t196 + pkin(1);
	t178 = -t201 * t193 - t202 * t194;
	t1 = [t181 * t189 + t186 * t190, t181 * t190 - t186 * t189, (-t196 * t209 + t208) * t197 + t198 * t204, t178 * t197 + t179 * t200 + 0; t180 * t189 - t190 * t185, t180 * t190 + t189 * t185, (t194 * t203 + t206) * t198 - t200 * t208, -t178 * t200 + t179 * t197 + 0; t184 * t189 - t190 * t210, t184 * t190 + t189 * t210, t194 * t195 + t196 * t205, t202 * t193 - t201 * t194 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end