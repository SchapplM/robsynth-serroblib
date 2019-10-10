% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPPRPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3,theta5]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:34
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPPRPR1_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR1_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRPR1_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRPR1_jacobiR_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:34:17
	% EndTime: 2019-10-09 23:34:17
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:34:17
	% EndTime: 2019-10-09 23:34:17
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
	% StartTime: 2019-10-09 23:34:17
	% EndTime: 2019-10-09 23:34:17
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (7->4), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->4)
	t12 = qJ(1) + pkin(9);
	t11 = cos(t12);
	t10 = sin(t12);
	t1 = [-t10, 0, 0, 0, 0, 0; t11, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t11, 0, 0, 0, 0, 0; -t10, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:34:17
	% EndTime: 2019-10-09 23:34:17
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->3), mult. (4->4), div. (0->0), fcn. (10->4), ass. (0->6)
	t12 = cos(pkin(10));
	t11 = sin(pkin(10));
	t10 = qJ(1) + pkin(9);
	t9 = cos(t10);
	t8 = sin(t10);
	t1 = [-t8 * t12, 0, 0, 0, 0, 0; t9 * t12, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t8 * t11, 0, 0, 0, 0, 0; -t9 * t11, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t9, 0, 0, 0, 0, 0; t8, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:34:17
	% EndTime: 2019-10-09 23:34:17
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (27->9), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->11)
	t18 = pkin(10) + qJ(4);
	t14 = sin(t18);
	t19 = qJ(1) + pkin(9);
	t15 = sin(t19);
	t23 = t15 * t14;
	t16 = cos(t18);
	t22 = t15 * t16;
	t17 = cos(t19);
	t21 = t17 * t14;
	t20 = t17 * t16;
	t1 = [-t22, 0, 0, -t21, 0, 0; t20, 0, 0, -t23, 0, 0; 0, 0, 0, t16, 0, 0; t23, 0, 0, -t20, 0, 0; -t21, 0, 0, -t22, 0, 0; 0, 0, 0, -t14, 0, 0; t17, 0, 0, 0, 0, 0; t15, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:34:17
	% EndTime: 2019-10-09 23:34:18
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (39->10), mult. (26->18), div. (0->0), fcn. (45->6), ass. (0->15)
	t59 = qJ(1) + pkin(9);
	t55 = sin(t59);
	t60 = sin(pkin(11));
	t67 = t55 * t60;
	t61 = cos(pkin(11));
	t66 = t55 * t61;
	t58 = pkin(10) + qJ(4);
	t56 = cos(t58);
	t65 = t56 * t60;
	t64 = t56 * t61;
	t57 = cos(t59);
	t63 = t57 * t60;
	t62 = t57 * t61;
	t54 = sin(t58);
	t1 = [-t55 * t64 + t63, 0, 0, -t54 * t62, 0, 0; t56 * t62 + t67, 0, 0, -t54 * t66, 0, 0; 0, 0, 0, t64, 0, 0; t55 * t65 + t62, 0, 0, t54 * t63, 0, 0; -t56 * t63 + t66, 0, 0, t54 * t67, 0, 0; 0, 0, 0, -t65, 0, 0; -t55 * t54, 0, 0, t57 * t56, 0, 0; t57 * t54, 0, 0, t55 * t56, 0, 0; 0, 0, 0, t54, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:34:18
	% EndTime: 2019-10-09 23:34:18
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (83->15), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->21)
	t80 = pkin(10) + qJ(4);
	t74 = sin(t80);
	t81 = qJ(1) + pkin(9);
	t75 = sin(t81);
	t88 = t75 * t74;
	t79 = pkin(11) + qJ(6);
	t76 = cos(t79);
	t87 = t75 * t76;
	t73 = sin(t79);
	t77 = cos(t80);
	t86 = t77 * t73;
	t85 = t77 * t76;
	t78 = cos(t81);
	t84 = t78 * t74;
	t83 = t78 * t76;
	t82 = t78 * t77;
	t72 = t75 * t73 + t76 * t82;
	t71 = -t73 * t82 + t87;
	t70 = t78 * t73 - t75 * t85;
	t69 = t75 * t86 + t83;
	t1 = [t70, 0, 0, -t74 * t83, 0, t71; t72, 0, 0, -t74 * t87, 0, -t69; 0, 0, 0, t85, 0, -t74 * t73; t69, 0, 0, t73 * t84, 0, -t72; t71, 0, 0, t73 * t88, 0, t70; 0, 0, 0, -t86, 0, -t74 * t76; -t88, 0, 0, t82, 0, 0; t84, 0, 0, t75 * t77, 0, 0; 0, 0, 0, t74, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end