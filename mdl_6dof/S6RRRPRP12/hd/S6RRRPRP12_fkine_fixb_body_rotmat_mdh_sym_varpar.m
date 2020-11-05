% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRPRP12 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:29
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRRPRP12_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP12_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRP12_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP12_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:29:37
	% EndTime: 2020-11-04 22:29:37
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:29:37
	% EndTime: 2020-11-04 22:29:37
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t98 = cos(qJ(1));
	t97 = sin(qJ(1));
	t1 = [t98, -t97, 0, 0; t97, t98, 0, 0; 0, 0, 1, pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:29:37
	% EndTime: 2020-11-04 22:29:37
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t102 = sin(qJ(1));
	t99 = sin(pkin(6));
	t110 = t102 * t99;
	t104 = cos(qJ(1));
	t109 = t104 * t99;
	t101 = sin(qJ(2));
	t108 = t102 * t101;
	t103 = cos(qJ(2));
	t107 = t102 * t103;
	t106 = t104 * t101;
	t105 = t104 * t103;
	t100 = cos(pkin(6));
	t1 = [-t100 * t108 + t105, -t100 * t107 - t106, t110, t104 * pkin(1) + pkin(8) * t110 + 0; t100 * t106 + t107, t100 * t105 - t108, -t109, t102 * pkin(1) - pkin(8) * t109 + 0; t99 * t101, t99 * t103, t100, t100 * pkin(8) + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:29:37
	% EndTime: 2020-11-04 22:29:37
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (29->23), mult. (60->37), div. (0->0), fcn. (81->8), ass. (0->20)
	t114 = sin(pkin(6));
	t119 = cos(qJ(3));
	t129 = t114 * t119;
	t116 = sin(qJ(3));
	t128 = t116 * t114;
	t117 = sin(qJ(2));
	t127 = t117 * t119;
	t118 = sin(qJ(1));
	t126 = t118 * t117;
	t120 = cos(qJ(2));
	t125 = t118 * t120;
	t121 = cos(qJ(1));
	t124 = t121 * t117;
	t123 = t121 * t120;
	t122 = pkin(2) * t117 - pkin(9) * t120;
	t115 = cos(pkin(6));
	t113 = t120 * pkin(2) + t117 * pkin(9) + pkin(1);
	t112 = -t115 * t124 - t125;
	t111 = t114 * pkin(8) - t122 * t115;
	t1 = [(-t115 * t127 + t128) * t118 + t119 * t123, (t115 * t126 - t123) * t116 + t118 * t129, t115 * t125 + t124, t111 * t118 + t113 * t121 + 0; -t112 * t119 - t121 * t128, t112 * t116 - t121 * t129, -t115 * t123 + t126, -t111 * t121 + t113 * t118 + 0; t114 * t127 + t115 * t116, t115 * t119 - t117 * t128, -t114 * t120, t115 * pkin(8) + t122 * t114 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:29:37
	% EndTime: 2020-11-04 22:29:37
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (45->27), mult. (78->41), div. (0->0), fcn. (99->8), ass. (0->22)
	t134 = sin(pkin(6));
	t139 = cos(qJ(3));
	t150 = t134 * t139;
	t136 = sin(qJ(3));
	t149 = t136 * t134;
	t137 = sin(qJ(2));
	t148 = t137 * t139;
	t138 = sin(qJ(1));
	t147 = t138 * t137;
	t140 = cos(qJ(2));
	t146 = t138 * t140;
	t141 = cos(qJ(1));
	t145 = t141 * t137;
	t144 = t141 * t140;
	t133 = t139 * pkin(3) + qJ(4) * t136 + pkin(2);
	t143 = t140 * pkin(9) - t133 * t137;
	t142 = t136 * pkin(3) - qJ(4) * t139 + pkin(8);
	t135 = cos(pkin(6));
	t132 = -t135 * t145 - t146;
	t131 = t137 * pkin(9) + t133 * t140 + pkin(1);
	t130 = t134 * t142 + t143 * t135;
	t1 = [t135 * t146 + t145, (t135 * t148 - t149) * t138 - t139 * t144, (-t135 * t147 + t144) * t136 - t138 * t150, t130 * t138 + t131 * t141 + 0; -t135 * t144 + t147, t132 * t139 + t141 * t149, -t132 * t136 + t141 * t150, -t130 * t141 + t131 * t138 + 0; -t134 * t140, -t134 * t148 - t135 * t136, -t135 * t139 + t137 * t149, -t143 * t134 + t142 * t135 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:29:37
	% EndTime: 2020-11-04 22:29:37
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (70->36), mult. (112->58), div. (0->0), fcn. (146->10), ass. (0->31)
	t158 = cos(pkin(6));
	t164 = cos(qJ(3));
	t180 = t158 * t164;
	t159 = sin(qJ(5));
	t165 = cos(qJ(2));
	t179 = t159 * t165;
	t157 = sin(pkin(6));
	t160 = sin(qJ(3));
	t178 = t160 * t157;
	t161 = sin(qJ(2));
	t177 = t160 * t161;
	t162 = sin(qJ(1));
	t176 = t162 * t165;
	t163 = cos(qJ(5));
	t175 = t163 * t165;
	t174 = t164 * t157;
	t166 = cos(qJ(1));
	t173 = t165 * t166;
	t172 = t166 * t161;
	t168 = pkin(3) + pkin(10);
	t156 = qJ(4) * t160 + t168 * t164 + pkin(2);
	t167 = pkin(4) + pkin(9);
	t171 = t156 * t161 - t167 * t165;
	t170 = qJ(4) * t164 - t168 * t160 - pkin(8);
	t153 = t158 * t177 + t174;
	t169 = -t153 * t159 + t158 * t175;
	t155 = t160 * t179 + t163 * t161;
	t154 = t157 * t177 - t180;
	t152 = t156 * t165 + t167 * t161 + pkin(1);
	t151 = -t157 * t170 - t171 * t158;
	t1 = [t166 * t155 + t169 * t162, (-t153 * t162 + t160 * t173) * t163 - (t158 * t176 + t172) * t159, (-t161 * t180 + t178) * t162 + t164 * t173, t151 * t162 + t152 * t166 + 0; t162 * t155 - t169 * t166, (t153 * t163 + t158 * t179) * t166 + t162 * (-t159 * t161 + t160 * t175), (t158 * t172 + t176) * t164 - t166 * t178, -t151 * t166 + t152 * t162 + 0; t154 * t159 - t157 * t175, t154 * t163 + t157 * t179, t158 * t160 + t161 * t174, t171 * t157 - t170 * t158 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:29:37
	% EndTime: 2020-11-04 22:29:37
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (96->40), mult. (138->62), div. (0->0), fcn. (172->10), ass. (0->32)
	t190 = cos(pkin(6));
	t196 = cos(qJ(3));
	t211 = t190 * t196;
	t191 = sin(qJ(5));
	t197 = cos(qJ(2));
	t210 = t191 * t197;
	t189 = sin(pkin(6));
	t192 = sin(qJ(3));
	t209 = t192 * t189;
	t193 = sin(qJ(2));
	t208 = t192 * t193;
	t194 = sin(qJ(1));
	t207 = t194 * t197;
	t195 = cos(qJ(5));
	t206 = t195 * t197;
	t205 = t196 * t189;
	t198 = cos(qJ(1));
	t204 = t197 * t198;
	t203 = t198 * t193;
	t188 = t191 * pkin(5) - qJ(6) * t195 + qJ(4);
	t199 = pkin(3) + pkin(10);
	t183 = t188 * t192 + t199 * t196 + pkin(2);
	t187 = t195 * pkin(5) + qJ(6) * t191 + pkin(4) + pkin(9);
	t202 = t183 * t193 - t187 * t197;
	t201 = t188 * t196 - t199 * t192 - pkin(8);
	t184 = t190 * t208 + t205;
	t200 = -t184 * t191 + t190 * t206;
	t186 = t192 * t210 + t195 * t193;
	t185 = t189 * t208 - t211;
	t182 = t183 * t197 + t187 * t193 + pkin(1);
	t181 = -t201 * t189 - t202 * t190;
	t1 = [t198 * t186 + t200 * t194, (-t193 * t211 + t209) * t194 + t196 * t204, (t184 * t194 - t192 * t204) * t195 + (t190 * t207 + t203) * t191, t181 * t194 + t182 * t198 + 0; t194 * t186 - t200 * t198, (t190 * t203 + t207) * t196 - t198 * t209, (-t184 * t195 - t190 * t210) * t198 - t194 * (-t191 * t193 + t192 * t206), -t181 * t198 + t182 * t194 + 0; t185 * t191 - t189 * t206, t190 * t192 + t193 * t205, -t185 * t195 - t189 * t210, t202 * t189 - t201 * t190 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end