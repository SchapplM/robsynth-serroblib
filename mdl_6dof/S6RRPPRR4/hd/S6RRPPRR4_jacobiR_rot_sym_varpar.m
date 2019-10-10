% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPPRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:41
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPPRR4_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR4_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_jacobiR_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:41:13
	% EndTime: 2019-10-10 09:41:13
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:41:13
	% EndTime: 2019-10-10 09:41:13
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
	% StartTime: 2019-10-10 09:41:13
	% EndTime: 2019-10-10 09:41:13
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
	% StartTime: 2019-10-10 09:41:13
	% EndTime: 2019-10-10 09:41:13
	% DurationCPUTime: 0.06s
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
	% StartTime: 2019-10-10 09:41:13
	% EndTime: 2019-10-10 09:41:13
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (26->8), mult. (74->18), div. (0->0), fcn. (112->8), ass. (0->17)
	t86 = cos(pkin(6));
	t83 = sin(pkin(11));
	t85 = cos(pkin(11));
	t87 = sin(qJ(2));
	t89 = cos(qJ(2));
	t91 = t89 * t83 + t87 * t85;
	t79 = t91 * t86;
	t80 = t87 * t83 - t89 * t85;
	t88 = sin(qJ(1));
	t90 = cos(qJ(1));
	t93 = t90 * t79 - t88 * t80;
	t92 = t88 * t79 + t90 * t80;
	t84 = sin(pkin(6));
	t78 = t80 * t86;
	t77 = -t88 * t78 + t90 * t91;
	t76 = -t90 * t78 - t88 * t91;
	t1 = [t90 * t84, 0, 0, 0, 0, 0; t88 * t84, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t93, t77, 0, 0, 0, 0; t92, -t76, 0, 0, 0, 0; 0, t80 * t84, 0, 0, 0, 0; t76, -t92, 0, 0, 0, 0; t77, t93, 0, 0, 0, 0; 0, t91 * t84, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:41:13
	% EndTime: 2019-10-10 09:41:13
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (63->16), mult. (182->36), div. (0->0), fcn. (266->10), ass. (0->27)
	t119 = sin(pkin(6));
	t124 = sin(qJ(1));
	t133 = t119 * t124;
	t127 = cos(qJ(1));
	t132 = t119 * t127;
	t118 = sin(pkin(11));
	t120 = cos(pkin(11));
	t123 = sin(qJ(2));
	t126 = cos(qJ(2));
	t114 = t123 * t118 - t126 * t120;
	t121 = cos(pkin(6));
	t128 = t114 * t121;
	t130 = t126 * t118 + t123 * t120;
	t107 = -t124 * t130 - t127 * t128;
	t122 = sin(qJ(5));
	t125 = cos(qJ(5));
	t131 = t107 * t122 + t125 * t132;
	t113 = t130 * t121;
	t106 = t127 * t113 - t124 * t114;
	t108 = -t124 * t113 - t127 * t114;
	t129 = -t107 * t125 + t122 * t132;
	t112 = t130 * t119;
	t111 = t114 * t119;
	t109 = t124 * t128 - t127 * t130;
	t105 = -t109 * t122 + t125 * t133;
	t104 = -t109 * t125 - t122 * t133;
	t1 = [t131, t108 * t122, 0, 0, t104, 0; t105, t106 * t122, 0, 0, t129, 0; 0, t112 * t122, 0, 0, t111 * t125 - t121 * t122, 0; -t129, t108 * t125, 0, 0, -t105, 0; t104, t106 * t125, 0, 0, t131, 0; 0, t112 * t125, 0, 0, -t111 * t122 - t121 * t125, 0; -t106, t109, 0, 0, 0, 0; t108, t107, 0, 0, 0, 0; 0, -t111, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:41:14
	% EndTime: 2019-10-10 09:41:14
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (157->33), mult. (425->65), div. (0->0), fcn. (606->12), ass. (0->37)
	t170 = sin(pkin(11));
	t172 = cos(pkin(11));
	t176 = sin(qJ(2));
	t180 = cos(qJ(2));
	t165 = t176 * t170 - t180 * t172;
	t173 = cos(pkin(6));
	t163 = t165 * t173;
	t177 = sin(qJ(1));
	t181 = cos(qJ(1));
	t182 = t180 * t170 + t176 * t172;
	t154 = -t181 * t163 - t177 * t182;
	t175 = sin(qJ(5));
	t179 = cos(qJ(5));
	t171 = sin(pkin(6));
	t189 = t171 * t181;
	t150 = t154 * t175 + t179 * t189;
	t174 = sin(qJ(6));
	t178 = cos(qJ(6));
	t164 = t182 * t173;
	t184 = t181 * t164 - t177 * t165;
	t194 = t150 * t174 + t184 * t178;
	t193 = t150 * t178 - t184 * t174;
	t190 = t171 * t177;
	t188 = t174 * t175;
	t187 = t175 * t178;
	t149 = -t154 * t179 + t175 * t189;
	t183 = -t177 * t164 - t181 * t165;
	t162 = t182 * t171;
	t161 = t165 * t171;
	t160 = t161 * t175 + t173 * t179;
	t159 = t161 * t179 - t173 * t175;
	t157 = t177 * t163 - t181 * t182;
	t148 = -t157 * t175 + t179 * t190;
	t147 = t157 * t179 + t175 * t190;
	t146 = t148 * t178 + t174 * t183;
	t145 = -t148 * t174 + t178 * t183;
	t1 = [t193, t157 * t174 + t183 * t187, 0, 0, -t147 * t178, t145; t146, t154 * t174 + t184 * t187, 0, 0, t149 * t178, t194; 0, -t161 * t174 + t162 * t187, 0, 0, t159 * t178, -t160 * t174 + t162 * t178; -t194, t157 * t178 - t183 * t188, 0, 0, t147 * t174, -t146; t145, t154 * t178 - t184 * t188, 0, 0, -t149 * t174, t193; 0, -t161 * t178 - t162 * t188, 0, 0, -t159 * t174, -t160 * t178 - t162 * t174; t149, -t183 * t179, 0, 0, t148, 0; t147, -t184 * t179, 0, 0, -t150, 0; 0, -t162 * t179, 0, 0, t160, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end