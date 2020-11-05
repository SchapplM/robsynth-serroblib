% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRPRR14 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:34
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRRPRR14_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR14_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR14_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR14_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:34:10
	% EndTime: 2020-11-04 22:34:10
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:34:10
	% EndTime: 2020-11-04 22:34:10
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t93 = cos(qJ(1));
	t92 = sin(qJ(1));
	t1 = [t93, -t92, 0, 0; t92, t93, 0, 0; 0, 0, 1, pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:34:10
	% EndTime: 2020-11-04 22:34:10
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t94 = sin(pkin(6));
	t97 = sin(qJ(1));
	t105 = t97 * t94;
	t96 = sin(qJ(2));
	t104 = t97 * t96;
	t98 = cos(qJ(2));
	t103 = t97 * t98;
	t99 = cos(qJ(1));
	t102 = t99 * t94;
	t101 = t99 * t96;
	t100 = t99 * t98;
	t95 = cos(pkin(6));
	t1 = [-t95 * t104 + t100, -t95 * t103 - t101, t105, t99 * pkin(1) + pkin(8) * t105 + 0; t95 * t101 + t103, t95 * t100 - t104, -t102, t97 * pkin(1) - pkin(8) * t102 + 0; t94 * t96, t94 * t98, t95, t95 * pkin(8) + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:34:10
	% EndTime: 2020-11-04 22:34:10
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (29->23), mult. (60->37), div. (0->0), fcn. (81->8), ass. (0->20)
	t109 = sin(pkin(6));
	t114 = cos(qJ(3));
	t124 = t109 * t114;
	t111 = sin(qJ(3));
	t123 = t111 * t109;
	t112 = sin(qJ(2));
	t122 = t112 * t114;
	t113 = sin(qJ(1));
	t121 = t113 * t112;
	t115 = cos(qJ(2));
	t120 = t113 * t115;
	t116 = cos(qJ(1));
	t119 = t116 * t112;
	t118 = t116 * t115;
	t117 = pkin(2) * t112 - pkin(9) * t115;
	t110 = cos(pkin(6));
	t108 = t115 * pkin(2) + t112 * pkin(9) + pkin(1);
	t107 = -t110 * t119 - t120;
	t106 = t109 * pkin(8) - t117 * t110;
	t1 = [(-t110 * t122 + t123) * t113 + t114 * t118, (t110 * t121 - t118) * t111 + t113 * t124, t110 * t120 + t119, t106 * t113 + t108 * t116 + 0; -t107 * t114 - t116 * t123, t107 * t111 - t116 * t124, -t110 * t118 + t121, -t106 * t116 + t108 * t113 + 0; t109 * t122 + t110 * t111, t110 * t114 - t112 * t123, -t109 * t115, t110 * pkin(8) + t117 * t109 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:34:10
	% EndTime: 2020-11-04 22:34:10
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (45->27), mult. (78->41), div. (0->0), fcn. (99->8), ass. (0->22)
	t129 = sin(pkin(6));
	t134 = cos(qJ(3));
	t145 = t129 * t134;
	t131 = sin(qJ(3));
	t144 = t131 * t129;
	t132 = sin(qJ(2));
	t143 = t132 * t134;
	t133 = sin(qJ(1));
	t142 = t133 * t132;
	t135 = cos(qJ(2));
	t141 = t133 * t135;
	t136 = cos(qJ(1));
	t140 = t136 * t132;
	t139 = t136 * t135;
	t128 = t134 * pkin(3) + qJ(4) * t131 + pkin(2);
	t138 = t135 * pkin(9) - t128 * t132;
	t137 = t131 * pkin(3) - qJ(4) * t134 + pkin(8);
	t130 = cos(pkin(6));
	t127 = -t130 * t140 - t141;
	t126 = t132 * pkin(9) + t128 * t135 + pkin(1);
	t125 = t129 * t137 + t138 * t130;
	t1 = [t130 * t141 + t140, (t130 * t143 - t144) * t133 - t134 * t139, (-t130 * t142 + t139) * t131 - t133 * t145, t125 * t133 + t126 * t136 + 0; -t130 * t139 + t142, t127 * t134 + t136 * t144, -t127 * t131 + t136 * t145, -t125 * t136 + t126 * t133 + 0; -t129 * t135, -t129 * t143 - t130 * t131, -t130 * t134 + t132 * t144, -t138 * t129 + t137 * t130 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:34:10
	% EndTime: 2020-11-04 22:34:10
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (70->36), mult. (112->58), div. (0->0), fcn. (146->10), ass. (0->31)
	t153 = cos(pkin(6));
	t159 = cos(qJ(3));
	t175 = t153 * t159;
	t154 = sin(qJ(5));
	t160 = cos(qJ(2));
	t174 = t154 * t160;
	t152 = sin(pkin(6));
	t155 = sin(qJ(3));
	t173 = t155 * t152;
	t156 = sin(qJ(2));
	t172 = t155 * t156;
	t157 = sin(qJ(1));
	t171 = t157 * t160;
	t158 = cos(qJ(5));
	t170 = t158 * t160;
	t169 = t159 * t152;
	t161 = cos(qJ(1));
	t168 = t160 * t161;
	t167 = t161 * t156;
	t163 = pkin(3) + pkin(10);
	t151 = qJ(4) * t155 + t163 * t159 + pkin(2);
	t162 = pkin(4) + pkin(9);
	t166 = t151 * t156 - t162 * t160;
	t165 = qJ(4) * t159 - t163 * t155 - pkin(8);
	t148 = t153 * t172 + t169;
	t164 = -t148 * t154 + t153 * t170;
	t150 = t155 * t174 + t158 * t156;
	t149 = t152 * t172 - t175;
	t147 = t151 * t160 + t162 * t156 + pkin(1);
	t146 = -t152 * t165 - t166 * t153;
	t1 = [t161 * t150 + t164 * t157, (-t148 * t157 + t155 * t168) * t158 - (t153 * t171 + t167) * t154, (-t156 * t175 + t173) * t157 + t159 * t168, t146 * t157 + t147 * t161 + 0; t157 * t150 - t164 * t161, (t148 * t158 + t153 * t174) * t161 + t157 * (-t154 * t156 + t155 * t170), (t153 * t167 + t171) * t159 - t161 * t173, -t146 * t161 + t147 * t157 + 0; t149 * t154 - t152 * t170, t149 * t158 + t152 * t174, t153 * t155 + t156 * t169, t166 * t152 - t165 * t153 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:34:10
	% EndTime: 2020-11-04 22:34:10
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (103->37), mult. (125->55), div. (0->0), fcn. (159->12), ass. (0->34)
	t191 = sin(pkin(6));
	t197 = cos(qJ(2));
	t208 = t191 * t197;
	t192 = cos(pkin(6));
	t196 = cos(qJ(3));
	t207 = t192 * t196;
	t193 = sin(qJ(3));
	t206 = t193 * t191;
	t194 = sin(qJ(2));
	t205 = t193 * t194;
	t195 = sin(qJ(1));
	t204 = t195 * t197;
	t203 = t196 * t191;
	t198 = cos(qJ(1));
	t202 = t197 * t198;
	t201 = t198 * t194;
	t186 = sin(qJ(5)) * pkin(5) + qJ(4);
	t189 = pkin(11) + pkin(3) + pkin(10);
	t180 = t186 * t193 + t189 * t196 + pkin(2);
	t185 = cos(qJ(5)) * pkin(5) + pkin(4) + pkin(9);
	t200 = t180 * t194 - t185 * t197;
	t199 = -t186 * t196 + t189 * t193 + pkin(8);
	t190 = qJ(5) + qJ(6);
	t188 = cos(t190);
	t187 = sin(t190);
	t184 = t192 * t204 + t201;
	t183 = t192 * t202 - t195 * t194;
	t182 = t191 * t205 - t207;
	t181 = t192 * t205 + t203;
	t179 = -t181 * t195 + t193 * t202;
	t178 = t181 * t198 + t193 * t204;
	t177 = t180 * t197 + t185 * t194 + pkin(1);
	t176 = t199 * t191 - t200 * t192;
	t1 = [t179 * t187 + t184 * t188, t179 * t188 - t184 * t187, (-t194 * t207 + t206) * t195 + t196 * t202, t176 * t195 + t177 * t198 + 0; t178 * t187 - t188 * t183, t178 * t188 + t187 * t183, (t192 * t201 + t204) * t196 - t198 * t206, -t176 * t198 + t177 * t195 + 0; t182 * t187 - t188 * t208, t182 * t188 + t187 * t208, t192 * t193 + t194 * t203, t200 * t191 + t199 * t192 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end