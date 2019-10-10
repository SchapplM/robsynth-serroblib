% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPPRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta4]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:53
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPPRR11_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR11_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR11_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR11_jacobiR_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:53:50
	% EndTime: 2019-10-10 09:53:50
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:53:50
	% EndTime: 2019-10-10 09:53:50
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
	% StartTime: 2019-10-10 09:53:50
	% EndTime: 2019-10-10 09:53:51
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
	% StartTime: 2019-10-10 09:53:51
	% EndTime: 2019-10-10 09:53:51
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->7), mult. (28->12), div. (0->0), fcn. (48->6), ass. (0->15)
	t62 = sin(qJ(2));
	t63 = sin(qJ(1));
	t69 = t63 * t62;
	t64 = cos(qJ(2));
	t68 = t63 * t64;
	t65 = cos(qJ(1));
	t67 = t65 * t62;
	t66 = t65 * t64;
	t61 = cos(pkin(6));
	t60 = sin(pkin(6));
	t59 = -t61 * t69 + t66;
	t58 = t61 * t68 + t67;
	t57 = t61 * t67 + t68;
	t56 = -t61 * t66 + t69;
	t1 = [t65 * t60, 0, 0, 0, 0, 0; t63 * t60, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t57, t58, 0, 0, 0, 0; -t59, t56, 0, 0, 0, 0; 0, -t60 * t64, 0, 0, 0, 0; -t56, t59, 0, 0, 0, 0; t58, t57, 0, 0, 0, 0; 0, t60 * t62, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:53:51
	% EndTime: 2019-10-10 09:53:51
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (16->10), mult. (57->26), div. (0->0), fcn. (88->8), ass. (0->20)
	t67 = sin(pkin(6));
	t70 = sin(qJ(2));
	t80 = t67 * t70;
	t71 = sin(qJ(1));
	t79 = t67 * t71;
	t73 = cos(qJ(1));
	t78 = t67 * t73;
	t77 = t71 * t70;
	t72 = cos(qJ(2));
	t76 = t71 * t72;
	t75 = t73 * t70;
	t74 = t73 * t72;
	t69 = cos(pkin(6));
	t68 = cos(pkin(11));
	t66 = sin(pkin(11));
	t65 = -t69 * t77 + t74;
	t64 = t69 * t76 + t75;
	t63 = t69 * t75 + t76;
	t62 = t69 * t74 - t77;
	t1 = [t62 * t66 + t68 * t78, t65 * t66, 0, 0, 0, 0; t64 * t66 + t68 * t79, t63 * t66, 0, 0, 0, 0; 0, t66 * t80, 0, 0, 0, 0; t62 * t68 - t66 * t78, t65 * t68, 0, 0, 0, 0; t64 * t68 - t66 * t79, t63 * t68, 0, 0, 0, 0; 0, t68 * t80, 0, 0, 0, 0; -t63, -t64, 0, 0, 0, 0; t65, t62, 0, 0, 0, 0; 0, t67 * t72, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:53:51
	% EndTime: 2019-10-10 09:53:51
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (52->16), mult. (87->30), div. (0->0), fcn. (134->8), ass. (0->26)
	t92 = sin(pkin(6));
	t94 = sin(qJ(2));
	t107 = t92 * t94;
	t95 = sin(qJ(1));
	t106 = t92 * t95;
	t96 = cos(qJ(2));
	t105 = t92 * t96;
	t97 = cos(qJ(1));
	t104 = t92 * t97;
	t103 = t95 * t94;
	t102 = t95 * t96;
	t101 = t97 * t94;
	t100 = t97 * t96;
	t93 = cos(pkin(6));
	t83 = -t93 * t100 + t103;
	t91 = pkin(11) + qJ(5);
	t89 = sin(t91);
	t90 = cos(t91);
	t99 = t90 * t104 - t83 * t89;
	t98 = t89 * t104 + t83 * t90;
	t86 = -t93 * t103 + t100;
	t85 = t93 * t102 + t101;
	t84 = t93 * t101 + t102;
	t82 = t90 * t106 + t85 * t89;
	t81 = -t89 * t106 + t85 * t90;
	t1 = [t99, t86 * t89, 0, 0, t81, 0; t82, t84 * t89, 0, 0, t98, 0; 0, t89 * t107, 0, 0, -t90 * t105 - t93 * t89, 0; -t98, t86 * t90, 0, 0, -t82, 0; t81, t84 * t90, 0, 0, t99, 0; 0, t90 * t107, 0, 0, t89 * t105 - t93 * t90, 0; -t84, -t85, 0, 0, 0, 0; t86, -t83, 0, 0, 0, 0; 0, t105, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:53:51
	% EndTime: 2019-10-10 09:53:51
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (128->32), mult. (219->63), div. (0->0), fcn. (320->10), ass. (0->37)
	t141 = cos(pkin(6));
	t146 = cos(qJ(2));
	t147 = cos(qJ(1));
	t148 = t147 * t146;
	t143 = sin(qJ(2));
	t144 = sin(qJ(1));
	t151 = t144 * t143;
	t131 = -t141 * t148 + t151;
	t139 = pkin(11) + qJ(5);
	t137 = sin(t139);
	t138 = cos(t139);
	t140 = sin(pkin(6));
	t154 = t140 * t147;
	t126 = -t131 * t137 + t138 * t154;
	t149 = t147 * t143;
	t150 = t144 * t146;
	t132 = t141 * t149 + t150;
	t142 = sin(qJ(6));
	t145 = cos(qJ(6));
	t162 = t126 * t142 + t132 * t145;
	t161 = t126 * t145 - t132 * t142;
	t158 = t137 * t142;
	t157 = t137 * t145;
	t156 = t140 * t144;
	t155 = t140 * t146;
	t153 = t142 * t143;
	t152 = t143 * t145;
	t125 = t131 * t138 + t137 * t154;
	t134 = -t141 * t151 + t148;
	t133 = t141 * t150 + t149;
	t130 = -t137 * t155 + t141 * t138;
	t129 = -t141 * t137 - t138 * t155;
	t124 = t133 * t137 + t138 * t156;
	t123 = -t133 * t138 + t137 * t156;
	t122 = t124 * t145 + t134 * t142;
	t121 = -t124 * t142 + t134 * t145;
	t1 = [t161, -t133 * t142 + t134 * t157, 0, 0, -t123 * t145, t121; t122, -t131 * t142 + t132 * t157, 0, 0, t125 * t145, t162; 0, (t137 * t152 + t142 * t146) * t140, 0, 0, t129 * t145, -t130 * t142 + t140 * t152; -t162, -t133 * t145 - t134 * t158, 0, 0, t123 * t142, -t122; t121, -t131 * t145 - t132 * t158, 0, 0, -t125 * t142, t161; 0, (-t137 * t153 + t145 * t146) * t140, 0, 0, -t129 * t142, -t130 * t145 - t140 * t153; t125, -t134 * t138, 0, 0, t124, 0; t123, -t132 * t138, 0, 0, -t126, 0; 0, -t140 * t143 * t138, 0, 0, t130, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end