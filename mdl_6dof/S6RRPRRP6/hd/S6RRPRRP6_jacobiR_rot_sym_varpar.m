% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPRRP6
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:37
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRRP6_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP6_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRP6_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP6_jacobiR_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:37:18
	% EndTime: 2019-10-10 10:37:18
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:37:18
	% EndTime: 2019-10-10 10:37:19
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
	% StartTime: 2019-10-10 10:37:19
	% EndTime: 2019-10-10 10:37:19
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
	% StartTime: 2019-10-10 10:37:19
	% EndTime: 2019-10-10 10:37:19
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (26->10), mult. (74->18), div. (0->0), fcn. (112->8), ass. (0->17)
	t70 = cos(pkin(6));
	t67 = sin(pkin(11));
	t69 = cos(pkin(11));
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
	% StartTime: 2019-10-10 10:37:19
	% EndTime: 2019-10-10 10:37:19
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (66->18), mult. (182->36), div. (0->0), fcn. (266->10), ass. (0->27)
	t115 = sin(pkin(6));
	t120 = sin(qJ(1));
	t128 = t115 * t120;
	t123 = cos(qJ(1));
	t127 = t115 * t123;
	t117 = cos(pkin(6));
	t114 = sin(pkin(11));
	t116 = cos(pkin(11));
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
	% StartTime: 2019-10-10 10:37:19
	% EndTime: 2019-10-10 10:37:19
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (154->31), mult. (425->65), div. (0->0), fcn. (606->12), ass. (0->37)
	t173 = sin(qJ(4));
	t177 = cos(qJ(4));
	t171 = cos(pkin(6));
	t168 = sin(pkin(11));
	t170 = cos(pkin(11));
	t174 = sin(qJ(2));
	t178 = cos(qJ(2));
	t181 = t178 * t168 + t174 * t170;
	t162 = t181 * t171;
	t163 = t174 * t168 - t178 * t170;
	t175 = sin(qJ(1));
	t179 = cos(qJ(1));
	t183 = t179 * t162 - t175 * t163;
	t169 = sin(pkin(6));
	t188 = t169 * t179;
	t148 = t173 * t188 - t177 * t183;
	t180 = t163 * t171;
	t152 = -t175 * t181 - t179 * t180;
	t172 = sin(qJ(5));
	t176 = cos(qJ(5));
	t193 = t148 * t172 - t152 * t176;
	t192 = t148 * t176 + t152 * t172;
	t189 = t169 * t175;
	t187 = t172 * t177;
	t185 = t176 * t177;
	t182 = -t175 * t162 - t179 * t163;
	t146 = -t173 * t183 - t177 * t188;
	t161 = t181 * t169;
	t160 = t163 * t169;
	t158 = t161 * t177 + t171 * t173;
	t157 = -t161 * t173 + t171 * t177;
	t155 = t175 * t180 - t179 * t181;
	t150 = t173 * t189 + t177 * t182;
	t149 = t173 * t182 - t177 * t189;
	t145 = t150 * t176 - t155 * t172;
	t144 = -t150 * t172 - t155 * t176;
	t1 = [t192, t155 * t185 + t172 * t182, 0, -t149 * t176, t144, 0; t145, t152 * t185 + t172 * t183, 0, t146 * t176, t193, 0; 0, -t160 * t185 + t161 * t172, 0, t157 * t176, -t158 * t172 + t160 * t176, 0; -t193, -t155 * t187 + t176 * t182, 0, t149 * t172, -t145, 0; t144, -t152 * t187 + t176 * t183, 0, -t146 * t172, t192, 0; 0, t160 * t187 + t161 * t176, 0, -t157 * t172, -t158 * t176 - t160 * t172, 0; t146, t155 * t173, 0, t150, 0, 0; t149, t152 * t173, 0, -t148, 0, 0; 0, -t160 * t173, 0, t158, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:37:19
	% EndTime: 2019-10-10 10:37:20
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (151->30), mult. (425->65), div. (0->0), fcn. (606->12), ass. (0->37)
	t202 = sin(qJ(4));
	t206 = cos(qJ(4));
	t200 = cos(pkin(6));
	t197 = sin(pkin(11));
	t199 = cos(pkin(11));
	t203 = sin(qJ(2));
	t207 = cos(qJ(2));
	t210 = t207 * t197 + t203 * t199;
	t191 = t210 * t200;
	t192 = t203 * t197 - t207 * t199;
	t204 = sin(qJ(1));
	t208 = cos(qJ(1));
	t212 = t208 * t191 - t204 * t192;
	t198 = sin(pkin(6));
	t217 = t198 * t208;
	t177 = t202 * t217 - t206 * t212;
	t209 = t192 * t200;
	t181 = -t204 * t210 - t208 * t209;
	t201 = sin(qJ(5));
	t205 = cos(qJ(5));
	t222 = t177 * t201 - t181 * t205;
	t221 = t177 * t205 + t181 * t201;
	t218 = t198 * t204;
	t216 = t201 * t206;
	t214 = t205 * t206;
	t211 = -t204 * t191 - t208 * t192;
	t175 = -t202 * t212 - t206 * t217;
	t190 = t210 * t198;
	t189 = t192 * t198;
	t187 = t190 * t206 + t200 * t202;
	t186 = -t190 * t202 + t200 * t206;
	t184 = t204 * t209 - t208 * t210;
	t179 = t202 * t218 + t206 * t211;
	t178 = t202 * t211 - t206 * t218;
	t174 = t179 * t205 - t184 * t201;
	t173 = t179 * t201 + t184 * t205;
	t1 = [t221, t184 * t214 + t201 * t211, 0, -t178 * t205, -t173, 0; t174, t181 * t214 + t201 * t212, 0, t175 * t205, t222, 0; 0, -t189 * t214 + t190 * t201, 0, t186 * t205, -t187 * t201 + t189 * t205, 0; t175, t184 * t202, 0, t179, 0, 0; t178, t181 * t202, 0, -t177, 0, 0; 0, -t189 * t202, 0, t187, 0, 0; t222, t184 * t216 - t205 * t211, 0, -t178 * t201, t174, 0; t173, t181 * t216 - t205 * t212, 0, t175 * t201, -t221, 0; 0, -t189 * t216 - t190 * t205, 0, t186 * t201, t187 * t205 + t189 * t201, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end