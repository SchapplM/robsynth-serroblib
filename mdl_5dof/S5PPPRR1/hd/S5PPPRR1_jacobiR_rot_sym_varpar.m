% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:17
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5PPPRR1_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR1_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPPRR1_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR1_jacobiR_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:17:16
	% EndTime: 2019-10-24 10:17:16
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:17:16
	% EndTime: 2019-10-24 10:17:16
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:17:16
	% EndTime: 2019-10-24 10:17:16
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:17:16
	% EndTime: 2019-10-24 10:17:16
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:17:16
	% EndTime: 2019-10-24 10:17:16
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (16->7), mult. (14->11), div. (0->0), fcn. (24->6), ass. (0->11)
	t26 = sin(pkin(7));
	t27 = cos(pkin(8));
	t31 = t26 * t27;
	t24 = pkin(9) + qJ(4);
	t22 = sin(t24);
	t28 = cos(pkin(7));
	t30 = t28 * t22;
	t23 = cos(t24);
	t29 = t28 * t23;
	t25 = sin(pkin(8));
	t1 = [0, 0, 0, t26 * t23 - t27 * t30, 0; 0, 0, 0, -t22 * t31 - t29, 0; 0, 0, 0, -t25 * t22, 0; 0, 0, 0, -t26 * t22 - t27 * t29, 0; 0, 0, 0, -t23 * t31 + t30, 0; 0, 0, 0, -t25 * t23, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:17:16
	% EndTime: 2019-10-24 10:17:16
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (44->14), mult. (57->30), div. (0->0), fcn. (88->8), ass. (0->19)
	t70 = sin(pkin(8));
	t74 = sin(qJ(5));
	t80 = t70 * t74;
	t75 = cos(qJ(5));
	t79 = t70 * t75;
	t71 = sin(pkin(7));
	t72 = cos(pkin(8));
	t78 = t71 * t72;
	t69 = pkin(9) + qJ(4);
	t67 = sin(t69);
	t73 = cos(pkin(7));
	t77 = t73 * t67;
	t68 = cos(t69);
	t76 = t73 * t68;
	t66 = t71 * t67 + t72 * t76;
	t65 = t71 * t68 - t72 * t77;
	t64 = t68 * t78 - t77;
	t63 = -t67 * t78 - t76;
	t1 = [0, 0, 0, t65 * t75, -t66 * t74 + t73 * t79; 0, 0, 0, t63 * t75, -t64 * t74 + t71 * t79; 0, 0, 0, -t67 * t79, -t68 * t80 - t72 * t75; 0, 0, 0, -t65 * t74, -t66 * t75 - t73 * t80; 0, 0, 0, -t63 * t74, -t64 * t75 - t71 * t80; 0, 0, 0, t67 * t80, -t68 * t79 + t72 * t74; 0, 0, 0, t66, 0; 0, 0, 0, t64, 0; 0, 0, 0, t70 * t68, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end