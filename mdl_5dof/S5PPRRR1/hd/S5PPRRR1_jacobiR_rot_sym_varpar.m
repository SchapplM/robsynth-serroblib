% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PPRRR1
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
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:20
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5PPRRR1_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR1_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPRRR1_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR1_jacobiR_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:20:29
	% EndTime: 2019-10-24 10:20:29
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:20:29
	% EndTime: 2019-10-24 10:20:29
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:20:29
	% EndTime: 2019-10-24 10:20:29
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:20:29
	% EndTime: 2019-10-24 10:20:29
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (11->6), mult. (4->4), div. (0->0), fcn. (10->4), ass. (0->6)
	t11 = cos(pkin(8));
	t10 = sin(pkin(8));
	t9 = pkin(9) + qJ(3);
	t8 = cos(t9);
	t7 = sin(t9);
	t1 = [0, 0, -t11 * t7, 0, 0; 0, 0, -t10 * t7, 0, 0; 0, 0, t8, 0, 0; 0, 0, -t11 * t8, 0, 0; 0, 0, -t10 * t8, 0, 0; 0, 0, -t7, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:20:29
	% EndTime: 2019-10-24 10:20:29
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (34->11), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->10)
	t16 = pkin(9) + qJ(3) + qJ(4);
	t14 = sin(t16);
	t17 = sin(pkin(8));
	t22 = t17 * t14;
	t15 = cos(t16);
	t21 = t17 * t15;
	t18 = cos(pkin(8));
	t20 = t18 * t14;
	t19 = t18 * t15;
	t1 = [0, 0, -t20, -t20, 0; 0, 0, -t22, -t22, 0; 0, 0, t15, t15, 0; 0, 0, -t19, -t19, 0; 0, 0, -t21, -t21, 0; 0, 0, -t14, -t14, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:20:29
	% EndTime: 2019-10-24 10:20:29
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (60->13), mult. (38->18), div. (0->0), fcn. (66->6), ass. (0->20)
	t77 = pkin(9) + qJ(3) + qJ(4);
	t76 = cos(t77);
	t80 = sin(qJ(5));
	t88 = t76 * t80;
	t78 = sin(pkin(8));
	t87 = t78 * t80;
	t81 = cos(qJ(5));
	t86 = t78 * t81;
	t79 = cos(pkin(8));
	t85 = t79 * t80;
	t84 = t79 * t81;
	t75 = sin(t77);
	t83 = t75 * t86;
	t82 = t75 * t84;
	t74 = t76 * t81;
	t73 = t79 * t76;
	t72 = t78 * t76;
	t71 = t75 * t85;
	t70 = t75 * t87;
	t1 = [0, 0, -t82, -t82, -t76 * t85 + t86; 0, 0, -t83, -t83, -t76 * t87 - t84; 0, 0, t74, t74, -t75 * t80; 0, 0, t71, t71, -t76 * t84 - t87; 0, 0, t70, t70, -t76 * t86 + t85; 0, 0, -t88, -t88, -t75 * t81; 0, 0, t73, t73, 0; 0, 0, t72, t72, 0; 0, 0, t75, t75, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end