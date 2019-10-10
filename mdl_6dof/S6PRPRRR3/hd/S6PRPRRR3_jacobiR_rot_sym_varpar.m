% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRPRRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:57
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRPRRR3_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR3_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRRR3_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR3_jacobiR_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:57:34
	% EndTime: 2019-10-09 21:57:34
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:57:34
	% EndTime: 2019-10-09 21:57:34
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:57:34
	% EndTime: 2019-10-09 21:57:34
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (5->5), mult. (14->12), div. (0->0), fcn. (24->6), ass. (0->9)
	t22 = cos(pkin(6));
	t23 = sin(qJ(2));
	t26 = t22 * t23;
	t24 = cos(qJ(2));
	t25 = t22 * t24;
	t21 = cos(pkin(11));
	t20 = sin(pkin(6));
	t19 = sin(pkin(11));
	t1 = [0, -t19 * t25 - t21 * t23, 0, 0, 0, 0; 0, -t19 * t23 + t21 * t25, 0, 0, 0, 0; 0, t20 * t24, 0, 0, 0, 0; 0, t19 * t26 - t21 * t24, 0, 0, 0, 0; 0, -t19 * t24 - t21 * t26, 0, 0, 0, 0; 0, -t20 * t23, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:57:35
	% EndTime: 2019-10-09 21:57:35
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->7), mult. (27->18), div. (0->0), fcn. (42->8), ass. (0->14)
	t55 = sin(pkin(6));
	t60 = cos(qJ(2));
	t63 = t55 * t60;
	t58 = cos(pkin(6));
	t59 = sin(qJ(2));
	t62 = t58 * t59;
	t61 = t58 * t60;
	t57 = cos(pkin(11));
	t56 = cos(pkin(12));
	t54 = sin(pkin(11));
	t53 = sin(pkin(12));
	t52 = -t54 * t61 - t57 * t59;
	t51 = -t54 * t59 + t57 * t61;
	t1 = [0, t52 * t56, 0, 0, 0, 0; 0, t51 * t56, 0, 0, 0, 0; 0, t56 * t63, 0, 0, 0, 0; 0, -t52 * t53, 0, 0, 0, 0; 0, -t51 * t53, 0, 0, 0, 0; 0, -t53 * t63, 0, 0, 0, 0; 0, -t54 * t62 + t57 * t60, 0, 0, 0, 0; 0, t54 * t60 + t57 * t62, 0, 0, 0, 0; 0, t55 * t59, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:57:35
	% EndTime: 2019-10-09 21:57:35
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (37->14), mult. (57->32), div. (0->0), fcn. (88->8), ass. (0->20)
	t72 = sin(pkin(11));
	t73 = sin(pkin(6));
	t83 = t72 * t73;
	t74 = cos(pkin(11));
	t82 = t73 * t74;
	t76 = sin(qJ(2));
	t81 = t73 * t76;
	t77 = cos(qJ(2));
	t80 = t73 * t77;
	t75 = cos(pkin(6));
	t79 = t75 * t76;
	t78 = t75 * t77;
	t71 = pkin(12) + qJ(4);
	t70 = cos(t71);
	t69 = sin(t71);
	t68 = -t72 * t79 + t74 * t77;
	t67 = -t72 * t78 - t74 * t76;
	t66 = t72 * t77 + t74 * t79;
	t65 = -t72 * t76 + t74 * t78;
	t1 = [0, t67 * t70, 0, -t68 * t69 + t70 * t83, 0, 0; 0, t65 * t70, 0, -t66 * t69 - t70 * t82, 0, 0; 0, t70 * t80, 0, -t69 * t81 + t75 * t70, 0, 0; 0, -t67 * t69, 0, -t68 * t70 - t69 * t83, 0, 0; 0, -t65 * t69, 0, -t66 * t70 + t69 * t82, 0, 0; 0, -t69 * t80, 0, -t75 * t69 - t70 * t81, 0, 0; 0, t68, 0, 0, 0, 0; 0, t66, 0, 0, 0, 0; 0, t81, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:57:35
	% EndTime: 2019-10-09 21:57:35
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (89->14), mult. (87->32), div. (0->0), fcn. (134->8), ass. (0->26)
	t94 = sin(pkin(11));
	t95 = sin(pkin(6));
	t105 = t94 * t95;
	t96 = cos(pkin(11));
	t104 = t95 * t96;
	t98 = sin(qJ(2));
	t103 = t95 * t98;
	t99 = cos(qJ(2));
	t102 = t95 * t99;
	t97 = cos(pkin(6));
	t101 = t97 * t98;
	t100 = t97 * t99;
	t93 = pkin(12) + qJ(4) + qJ(5);
	t92 = cos(t93);
	t91 = sin(t93);
	t90 = -t94 * t101 + t96 * t99;
	t89 = -t94 * t100 - t96 * t98;
	t88 = t96 * t101 + t94 * t99;
	t87 = t96 * t100 - t94 * t98;
	t86 = -t92 * t103 - t97 * t91;
	t85 = -t91 * t103 + t97 * t92;
	t84 = -t91 * t105 - t90 * t92;
	t83 = t92 * t105 - t90 * t91;
	t82 = t91 * t104 - t88 * t92;
	t81 = -t92 * t104 - t88 * t91;
	t1 = [0, t89 * t92, 0, t83, t83, 0; 0, t87 * t92, 0, t81, t81, 0; 0, t92 * t102, 0, t85, t85, 0; 0, -t89 * t91, 0, t84, t84, 0; 0, -t87 * t91, 0, t82, t82, 0; 0, -t91 * t102, 0, t86, t86, 0; 0, t90, 0, 0, 0, 0; 0, t88, 0, 0, 0, 0; 0, t103, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:57:35
	% EndTime: 2019-10-09 21:57:35
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (186->31), mult. (214->65), div. (0->0), fcn. (313->10), ass. (0->37)
	t146 = sin(pkin(11));
	t148 = cos(pkin(11));
	t153 = cos(qJ(2));
	t149 = cos(pkin(6));
	t151 = sin(qJ(2));
	t157 = t149 * t151;
	t139 = t146 * t153 + t148 * t157;
	t145 = pkin(12) + qJ(4) + qJ(5);
	t143 = sin(t145);
	t144 = cos(t145);
	t147 = sin(pkin(6));
	t159 = t147 * t148;
	t131 = -t139 * t143 - t144 * t159;
	t150 = sin(qJ(6));
	t165 = t131 * t150;
	t141 = -t146 * t157 + t148 * t153;
	t160 = t146 * t147;
	t133 = -t141 * t143 + t144 * t160;
	t164 = t133 * t150;
	t158 = t147 * t151;
	t136 = -t143 * t158 + t149 * t144;
	t163 = t136 * t150;
	t162 = t144 * t150;
	t152 = cos(qJ(6));
	t161 = t144 * t152;
	t156 = t149 * t153;
	t155 = t150 * t153;
	t154 = t152 * t153;
	t140 = t146 * t156 + t148 * t151;
	t138 = t146 * t151 - t148 * t156;
	t137 = t149 * t143 + t144 * t158;
	t135 = t136 * t152;
	t134 = t141 * t144 + t143 * t160;
	t132 = t139 * t144 - t143 * t159;
	t130 = t133 * t152;
	t129 = t131 * t152;
	t1 = [0, -t140 * t161 + t141 * t150, 0, t130, t130, -t134 * t150 + t140 * t152; 0, -t138 * t161 + t139 * t150, 0, t129, t129, -t132 * t150 + t138 * t152; 0, (t144 * t154 + t150 * t151) * t147, 0, t135, t135, -t137 * t150 - t147 * t154; 0, t140 * t162 + t141 * t152, 0, -t164, -t164, -t134 * t152 - t140 * t150; 0, t138 * t162 + t139 * t152, 0, -t165, -t165, -t132 * t152 - t138 * t150; 0, (-t144 * t155 + t151 * t152) * t147, 0, -t163, -t163, -t137 * t152 + t147 * t155; 0, -t140 * t143, 0, t134, t134, 0; 0, -t138 * t143, 0, t132, t132, 0; 0, t147 * t153 * t143, 0, t137, t137, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end