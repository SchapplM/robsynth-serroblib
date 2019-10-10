% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:55
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRRR4_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR4_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRR4_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR4_jacobiR_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:55:48
	% EndTime: 2019-10-10 10:55:48
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:55:48
	% EndTime: 2019-10-10 10:55:48
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t8, 0, 0, 0, 0, 0; t9, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t9, 0, 0, 0, 0, 0; -t8, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:55:48
	% EndTime: 2019-10-10 10:55:48
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->7), mult. (28->12), div. (0->0), fcn. (48->6), ass. (0->15)
	t50 = sin(qJ(2));
	t51 = sin(qJ(1));
	t57 = t51 * t50;
	t52 = cos(qJ(2));
	t56 = t51 * t52;
	t53 = cos(qJ(1));
	t55 = t53 * t50;
	t54 = t53 * t52;
	t49 = cos(pkin(6));
	t48 = sin(pkin(6));
	t47 = -t49 * t57 + t54;
	t46 = -t49 * t56 - t55;
	t45 = -t49 * t55 - t56;
	t44 = -t49 * t54 + t57;
	t1 = [t45, t46, 0, 0, 0, 0; t47, -t44, 0, 0, 0, 0; 0, t48 * t52, 0, 0, 0, 0; t44, -t47, 0, 0, 0, 0; t46, t45, 0, 0, 0, 0; 0, -t48 * t50, 0, 0, 0, 0; t53 * t48, 0, 0, 0, 0, 0; t51 * t48, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:55:48
	% EndTime: 2019-10-10 10:55:48
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (26->10), mult. (74->18), div. (0->0), fcn. (112->8), ass. (0->17)
	t70 = cos(pkin(6));
	t67 = sin(pkin(12));
	t69 = cos(pkin(12));
	t71 = sin(qJ(2));
	t73 = cos(qJ(2));
	t75 = t73 * t67 + t71 * t69;
	t63 = t75 * t70;
	t64 = t71 * t67 - t73 * t69;
	t72 = sin(qJ(1));
	t74 = cos(qJ(1));
	t77 = -t74 * t63 + t72 * t64;
	t76 = t72 * t63 + t74 * t64;
	t68 = sin(pkin(6));
	t62 = t64 * t70;
	t61 = t72 * t62 - t74 * t75;
	t60 = -t74 * t62 - t72 * t75;
	t1 = [t77, t61, 0, 0, 0, 0; -t76, t60, 0, 0, 0, 0; 0, -t64 * t68, 0, 0, 0, 0; -t60, t76, 0, 0, 0, 0; t61, t77, 0, 0, 0, 0; 0, -t75 * t68, 0, 0, 0, 0; t74 * t68, 0, 0, 0, 0, 0; t72 * t68, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:55:48
	% EndTime: 2019-10-10 10:55:48
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (66->18), mult. (182->36), div. (0->0), fcn. (266->10), ass. (0->27)
	t115 = sin(pkin(6));
	t120 = sin(qJ(1));
	t128 = t115 * t120;
	t123 = cos(qJ(1));
	t127 = t115 * t123;
	t117 = cos(pkin(6));
	t114 = sin(pkin(12));
	t116 = cos(pkin(12));
	t119 = sin(qJ(2));
	t122 = cos(qJ(2));
	t125 = t122 * t114 + t119 * t116;
	t109 = t125 * t117;
	t110 = t119 * t114 - t122 * t116;
	t103 = t123 * t109 - t120 * t110;
	t118 = sin(qJ(4));
	t121 = cos(qJ(4));
	t126 = -t103 * t121 + t118 * t127;
	t105 = -t120 * t109 - t123 * t110;
	t124 = t103 * t118 + t121 * t127;
	t108 = t110 * t117;
	t107 = t125 * t115;
	t106 = t110 * t115;
	t104 = t120 * t108 - t123 * t125;
	t102 = -t123 * t108 - t120 * t125;
	t101 = t105 * t121 + t118 * t128;
	t100 = -t105 * t118 + t121 * t128;
	t1 = [t126, t104 * t121, 0, t100, 0, 0; t101, t102 * t121, 0, -t124, 0, 0; 0, -t106 * t121, 0, -t107 * t118 + t117 * t121, 0, 0; t124, -t104 * t118, 0, -t101, 0, 0; t100, -t102 * t118, 0, t126, 0, 0; 0, t106 * t118, 0, -t107 * t121 - t117 * t118, 0, 0; t102, t105, 0, 0, 0, 0; -t104, t103, 0, 0, 0, 0; 0, t107, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:55:48
	% EndTime: 2019-10-10 10:55:49
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (124->21), mult. (238->36), div. (0->0), fcn. (348->10), ass. (0->30)
	t140 = sin(pkin(6));
	t144 = sin(qJ(1));
	t150 = t140 * t144;
	t146 = cos(qJ(1));
	t149 = t140 * t146;
	t142 = cos(pkin(6));
	t139 = sin(pkin(12));
	t141 = cos(pkin(12));
	t143 = sin(qJ(2));
	t145 = cos(qJ(2));
	t148 = t145 * t139 + t143 * t141;
	t131 = t148 * t142;
	t132 = t143 * t139 - t145 * t141;
	t123 = t146 * t131 - t144 * t132;
	t138 = qJ(4) + qJ(5);
	t136 = sin(t138);
	t137 = cos(t138);
	t119 = -t123 * t137 + t136 * t149;
	t125 = -t144 * t131 - t146 * t132;
	t147 = t123 * t136 + t137 * t149;
	t130 = t132 * t142;
	t129 = t148 * t140;
	t128 = t132 * t140;
	t127 = -t129 * t137 - t142 * t136;
	t126 = -t129 * t136 + t142 * t137;
	t124 = t144 * t130 - t146 * t148;
	t122 = -t146 * t130 - t144 * t148;
	t121 = t125 * t137 + t136 * t150;
	t120 = -t125 * t136 + t137 * t150;
	t1 = [t119, t124 * t137, 0, t120, t120, 0; t121, t122 * t137, 0, -t147, -t147, 0; 0, -t128 * t137, 0, t126, t126, 0; t147, -t124 * t136, 0, -t121, -t121, 0; t120, -t122 * t136, 0, t119, t119, 0; 0, t128 * t136, 0, t127, t127, 0; t122, t125, 0, 0, 0, 0; -t124, t123, 0, 0, 0, 0; 0, t129, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:55:49
	% EndTime: 2019-10-10 10:55:49
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (256->36), mult. (515->65), div. (0->0), fcn. (735->12), ass. (0->44)
	t199 = qJ(4) + qJ(5);
	t197 = sin(t199);
	t198 = cos(t199);
	t203 = cos(pkin(6));
	t200 = sin(pkin(12));
	t202 = cos(pkin(12));
	t205 = sin(qJ(2));
	t208 = cos(qJ(2));
	t211 = t208 * t200 + t205 * t202;
	t191 = t211 * t203;
	t192 = t205 * t200 - t208 * t202;
	t206 = sin(qJ(1));
	t209 = cos(qJ(1));
	t213 = t209 * t191 - t206 * t192;
	t201 = sin(pkin(6));
	t216 = t201 * t209;
	t176 = t197 * t216 - t198 * t213;
	t210 = t192 * t203;
	t180 = -t206 * t211 - t209 * t210;
	t204 = sin(qJ(6));
	t207 = cos(qJ(6));
	t226 = t176 * t204 - t180 * t207;
	t225 = t176 * t207 + t180 * t204;
	t174 = -t197 * t213 - t198 * t216;
	t224 = t174 * t204;
	t212 = -t206 * t191 - t209 * t192;
	t217 = t201 * t206;
	t177 = t197 * t212 - t198 * t217;
	t223 = t177 * t204;
	t190 = t211 * t201;
	t186 = -t190 * t197 + t203 * t198;
	t220 = t186 * t204;
	t219 = t198 * t204;
	t218 = t198 * t207;
	t189 = t192 * t201;
	t187 = t190 * t198 + t203 * t197;
	t185 = t186 * t207;
	t183 = t206 * t210 - t209 * t211;
	t178 = t197 * t217 + t198 * t212;
	t173 = t177 * t207;
	t172 = t174 * t207;
	t171 = t178 * t207 - t183 * t204;
	t170 = -t178 * t204 - t183 * t207;
	t1 = [t225, t183 * t218 + t204 * t212, 0, -t173, -t173, t170; t171, t180 * t218 + t204 * t213, 0, t172, t172, t226; 0, -t189 * t218 + t190 * t204, 0, t185, t185, -t187 * t204 + t189 * t207; -t226, -t183 * t219 + t207 * t212, 0, t223, t223, -t171; t170, -t180 * t219 + t207 * t213, 0, -t224, -t224, t225; 0, t189 * t219 + t190 * t207, 0, -t220, -t220, -t187 * t207 - t189 * t204; t174, t183 * t197, 0, t178, t178, 0; t177, t180 * t197, 0, -t176, -t176, 0; 0, -t189 * t197, 0, t187, t187, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end