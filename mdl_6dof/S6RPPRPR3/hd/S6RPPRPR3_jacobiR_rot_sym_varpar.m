% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPPRPR3
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta5]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:37
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPPRPR3_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR3_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRPR3_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR3_jacobiR_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:37:42
	% EndTime: 2019-10-09 23:37:42
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:37:42
	% EndTime: 2019-10-09 23:37:42
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
	% StartTime: 2019-10-09 23:37:42
	% EndTime: 2019-10-09 23:37:42
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
	% StartTime: 2019-10-09 23:37:43
	% EndTime: 2019-10-09 23:37:43
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (5->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->4)
	t8 = qJ(1) + pkin(9);
	t7 = cos(t8);
	t6 = sin(t8);
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t6, 0, 0, 0, 0, 0; -t7, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t7, 0, 0, 0, 0, 0; t6, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:37:42
	% EndTime: 2019-10-09 23:37:43
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (15->6), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->10)
	t15 = qJ(1) + pkin(9);
	t13 = sin(t15);
	t16 = sin(qJ(4));
	t19 = t13 * t16;
	t14 = cos(t15);
	t17 = cos(qJ(4));
	t18 = t14 * t17;
	t12 = t14 * t16;
	t11 = t13 * t17;
	t1 = [t12, 0, 0, t11, 0, 0; t19, 0, 0, -t18, 0, 0; 0, 0, 0, -t16, 0, 0; t18, 0, 0, -t19, 0, 0; t11, 0, 0, t12, 0, 0; 0, 0, 0, -t17, 0, 0; -t13, 0, 0, 0, 0, 0; t14, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:37:42
	% EndTime: 2019-10-09 23:37:42
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (25->7), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->11)
	t21 = qJ(4) + pkin(10);
	t17 = sin(t21);
	t22 = qJ(1) + pkin(9);
	t18 = sin(t22);
	t24 = t18 * t17;
	t19 = cos(t21);
	t20 = cos(t22);
	t23 = t20 * t19;
	t16 = t20 * t17;
	t15 = t18 * t19;
	t1 = [t16, 0, 0, t15, 0, 0; t24, 0, 0, -t23, 0, 0; 0, 0, 0, -t17, 0, 0; t23, 0, 0, -t24, 0, 0; t15, 0, 0, t16, 0, 0; 0, 0, 0, -t19, 0, 0; -t18, 0, 0, 0, 0, 0; t20, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:37:43
	% EndTime: 2019-10-09 23:37:43
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (61->16), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->19)
	t80 = qJ(1) + pkin(9);
	t76 = sin(t80);
	t81 = sin(qJ(6));
	t88 = t76 * t81;
	t82 = cos(qJ(6));
	t87 = t76 * t82;
	t79 = qJ(4) + pkin(10);
	t77 = cos(t79);
	t86 = t77 * t81;
	t85 = t77 * t82;
	t78 = cos(t80);
	t84 = t78 * t81;
	t83 = t78 * t82;
	t75 = sin(t79);
	t74 = t75 * t83 - t88;
	t73 = t75 * t84 + t87;
	t72 = t75 * t87 + t84;
	t71 = -t75 * t88 + t83;
	t1 = [t74, 0, 0, t76 * t85, 0, t71; t72, 0, 0, -t77 * t83, 0, t73; 0, 0, 0, -t75 * t82, 0, -t86; -t73, 0, 0, -t76 * t86, 0, -t72; t71, 0, 0, t77 * t84, 0, t74; 0, 0, 0, t75 * t81, 0, -t85; -t78 * t77, 0, 0, t76 * t75, 0, 0; -t76 * t77, 0, 0, -t78 * t75, 0, 0; 0, 0, 0, t77, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end