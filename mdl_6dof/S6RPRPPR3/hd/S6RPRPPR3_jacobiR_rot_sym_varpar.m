% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRPPR3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:18
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRPPR3_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR3_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPPR3_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR3_jacobiR_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:18:44
	% EndTime: 2019-10-10 00:18:44
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:18:43
	% EndTime: 2019-10-10 00:18:43
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
	% StartTime: 2019-10-10 00:18:43
	% EndTime: 2019-10-10 00:18:44
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
	% StartTime: 2019-10-10 00:18:43
	% EndTime: 2019-10-10 00:18:44
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
	% StartTime: 2019-10-10 00:18:43
	% EndTime: 2019-10-10 00:18:43
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (14->5), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->10)
	t43 = qJ(1) + pkin(9);
	t41 = sin(t43);
	t44 = sin(qJ(3));
	t48 = t41 * t44;
	t45 = cos(qJ(3));
	t47 = t41 * t45;
	t42 = cos(t43);
	t46 = t42 * t44;
	t40 = t42 * t45;
	t1 = [-t47, 0, -t46, 0, 0, 0; t40, 0, -t48, 0, 0, 0; 0, 0, t45, 0, 0, 0; t42, 0, 0, 0, 0, 0; t41, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t48, 0, t40, 0, 0, 0; t46, 0, t47, 0, 0, 0; 0, 0, t44, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:18:43
	% EndTime: 2019-10-10 00:18:43
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (15->6), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->10)
	t17 = qJ(1) + pkin(9);
	t15 = sin(t17);
	t18 = sin(qJ(3));
	t21 = t15 * t18;
	t16 = cos(t17);
	t19 = cos(qJ(3));
	t20 = t16 * t19;
	t14 = t16 * t18;
	t13 = t15 * t19;
	t1 = [-t21, 0, t20, 0, 0, 0; t14, 0, t13, 0, 0, 0; 0, 0, t18, 0, 0, 0; t13, 0, t14, 0, 0, 0; -t20, 0, t21, 0, 0, 0; 0, 0, -t19, 0, 0, 0; -t16, 0, 0, 0, 0, 0; -t15, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:18:44
	% EndTime: 2019-10-10 00:18:44
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (38->13), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->16)
	t68 = sin(qJ(6));
	t69 = sin(qJ(3));
	t75 = t69 * t68;
	t70 = cos(qJ(6));
	t74 = t69 * t70;
	t71 = cos(qJ(3));
	t73 = t71 * t68;
	t72 = t71 * t70;
	t67 = qJ(1) + pkin(9);
	t66 = cos(t67);
	t65 = sin(t67);
	t64 = -t65 * t68 + t66 * t74;
	t63 = -t65 * t70 - t66 * t75;
	t62 = -t65 * t74 - t66 * t68;
	t61 = t65 * t75 - t66 * t70;
	t1 = [t62, 0, t66 * t72, 0, 0, t63; t64, 0, t65 * t72, 0, 0, -t61; 0, 0, t74, 0, 0, t73; t61, 0, -t66 * t73, 0, 0, -t64; t63, 0, -t65 * t73, 0, 0, t62; 0, 0, -t75, 0, 0, t72; -t65 * t71, 0, -t66 * t69, 0, 0, 0; t66 * t71, 0, -t65 * t69, 0, 0, 0; 0, 0, t71, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end