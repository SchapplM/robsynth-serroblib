% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRPPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4,theta5]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:15
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRPPR1_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR1_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPPR1_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPPR1_jacobiR_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:15:13
	% EndTime: 2019-10-10 00:15:13
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:15:13
	% EndTime: 2019-10-10 00:15:13
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
	% StartTime: 2019-10-10 00:15:13
	% EndTime: 2019-10-10 00:15:13
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
	% StartTime: 2019-10-10 00:15:13
	% EndTime: 2019-10-10 00:15:13
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (17->8), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->10)
	t13 = qJ(1) + pkin(9);
	t11 = sin(t13);
	t14 = sin(qJ(3));
	t19 = t11 * t14;
	t15 = cos(qJ(3));
	t18 = t11 * t15;
	t12 = cos(t13);
	t17 = t12 * t14;
	t16 = t12 * t15;
	t1 = [-t18, 0, -t17, 0, 0, 0; t16, 0, -t19, 0, 0, 0; 0, 0, t15, 0, 0, 0; t19, 0, -t16, 0, 0, 0; -t17, 0, -t18, 0, 0, 0; 0, 0, -t14, 0, 0, 0; t12, 0, 0, 0, 0, 0; t11, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:15:13
	% EndTime: 2019-10-10 00:15:14
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (27->9), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->11)
	t21 = qJ(3) + pkin(10);
	t17 = sin(t21);
	t22 = qJ(1) + pkin(9);
	t18 = sin(t22);
	t26 = t18 * t17;
	t19 = cos(t21);
	t25 = t18 * t19;
	t20 = cos(t22);
	t24 = t20 * t17;
	t23 = t20 * t19;
	t1 = [-t25, 0, -t24, 0, 0, 0; t23, 0, -t26, 0, 0, 0; 0, 0, t19, 0, 0, 0; t26, 0, -t23, 0, 0, 0; -t24, 0, -t25, 0, 0, 0; 0, 0, -t17, 0, 0, 0; t20, 0, 0, 0, 0, 0; t18, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:15:14
	% EndTime: 2019-10-10 00:15:14
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (39->10), mult. (26->18), div. (0->0), fcn. (45->6), ass. (0->15)
	t62 = qJ(1) + pkin(9);
	t58 = sin(t62);
	t63 = sin(pkin(11));
	t70 = t58 * t63;
	t64 = cos(pkin(11));
	t69 = t58 * t64;
	t61 = qJ(3) + pkin(10);
	t59 = cos(t61);
	t68 = t59 * t63;
	t67 = t59 * t64;
	t60 = cos(t62);
	t66 = t60 * t63;
	t65 = t60 * t64;
	t57 = sin(t61);
	t1 = [-t58 * t67 + t66, 0, -t57 * t65, 0, 0, 0; t59 * t65 + t70, 0, -t57 * t69, 0, 0, 0; 0, 0, t67, 0, 0, 0; t58 * t68 + t65, 0, t57 * t66, 0, 0, 0; -t59 * t66 + t69, 0, t57 * t70, 0, 0, 0; 0, 0, -t68, 0, 0, 0; -t58 * t57, 0, t60 * t59, 0, 0, 0; t60 * t57, 0, t58 * t59, 0, 0, 0; 0, 0, t57, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:15:13
	% EndTime: 2019-10-10 00:15:13
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (83->15), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->21)
	t83 = qJ(3) + pkin(10);
	t77 = sin(t83);
	t84 = qJ(1) + pkin(9);
	t78 = sin(t84);
	t91 = t78 * t77;
	t82 = pkin(11) + qJ(6);
	t79 = cos(t82);
	t90 = t78 * t79;
	t76 = sin(t82);
	t80 = cos(t83);
	t89 = t80 * t76;
	t88 = t80 * t79;
	t81 = cos(t84);
	t87 = t81 * t77;
	t86 = t81 * t79;
	t85 = t81 * t80;
	t75 = t78 * t76 + t79 * t85;
	t74 = -t76 * t85 + t90;
	t73 = t81 * t76 - t78 * t88;
	t72 = t78 * t89 + t86;
	t1 = [t73, 0, -t77 * t86, 0, 0, t74; t75, 0, -t77 * t90, 0, 0, -t72; 0, 0, t88, 0, 0, -t77 * t76; t72, 0, t76 * t87, 0, 0, -t75; t74, 0, t76 * t91, 0, 0, t73; 0, 0, -t89, 0, 0, -t77 * t79; -t91, 0, t85, 0, 0, 0; t87, 0, t78 * t80, 0, 0, 0; 0, 0, t77, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end