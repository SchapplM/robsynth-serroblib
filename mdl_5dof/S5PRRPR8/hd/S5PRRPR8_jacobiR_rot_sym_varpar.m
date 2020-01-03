% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PRRPR8
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
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5PRRPR8_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR8_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRPR8_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR8_jacobiR_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:43:02
	% EndTime: 2019-12-31 17:43:02
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:43:02
	% EndTime: 2019-12-31 17:43:02
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:43:02
	% EndTime: 2019-12-31 17:43:02
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (5->5), mult. (4->4), div. (0->0), fcn. (10->4), ass. (0->5)
	t9 = cos(qJ(2));
	t8 = sin(qJ(2));
	t7 = cos(pkin(8));
	t6 = sin(pkin(8));
	t1 = [0, -t7 * t8, 0, 0, 0; 0, -t6 * t8, 0, 0, 0; 0, t9, 0, 0, 0; 0, -t7 * t9, 0, 0, 0; 0, -t6 * t9, 0, 0, 0; 0, -t8, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:43:02
	% EndTime: 2019-12-31 17:43:02
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (22->11), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->10)
	t15 = qJ(2) + qJ(3);
	t13 = sin(t15);
	t16 = sin(pkin(8));
	t21 = t16 * t13;
	t14 = cos(t15);
	t20 = t16 * t14;
	t17 = cos(pkin(8));
	t19 = t17 * t13;
	t18 = t17 * t14;
	t1 = [0, -t19, -t19, 0, 0; 0, -t21, -t21, 0, 0; 0, t14, t14, 0, 0; 0, -t18, -t18, 0, 0; 0, -t20, -t20, 0, 0; 0, -t13, -t13, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:43:02
	% EndTime: 2019-12-31 17:43:02
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (34->11), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->10)
	t20 = qJ(2) + qJ(3) + pkin(9);
	t18 = sin(t20);
	t21 = sin(pkin(8));
	t26 = t21 * t18;
	t19 = cos(t20);
	t25 = t21 * t19;
	t22 = cos(pkin(8));
	t24 = t22 * t18;
	t23 = t22 * t19;
	t1 = [0, -t24, -t24, 0, 0; 0, -t26, -t26, 0, 0; 0, t19, t19, 0, 0; 0, -t23, -t23, 0, 0; 0, -t25, -t25, 0, 0; 0, -t18, -t18, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:43:02
	% EndTime: 2019-12-31 17:43:02
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (60->13), mult. (38->18), div. (0->0), fcn. (66->6), ass. (0->20)
	t81 = qJ(2) + qJ(3) + pkin(9);
	t80 = cos(t81);
	t84 = sin(qJ(5));
	t92 = t80 * t84;
	t82 = sin(pkin(8));
	t91 = t82 * t84;
	t85 = cos(qJ(5));
	t90 = t82 * t85;
	t83 = cos(pkin(8));
	t89 = t83 * t84;
	t88 = t83 * t85;
	t79 = sin(t81);
	t87 = t79 * t90;
	t86 = t79 * t88;
	t78 = t80 * t85;
	t77 = t83 * t80;
	t76 = t82 * t80;
	t75 = t79 * t89;
	t74 = t79 * t91;
	t1 = [0, -t86, -t86, 0, -t80 * t89 + t90; 0, -t87, -t87, 0, -t80 * t91 - t88; 0, t78, t78, 0, -t79 * t84; 0, t75, t75, 0, -t80 * t88 - t91; 0, t74, t74, 0, -t80 * t90 + t89; 0, -t92, -t92, 0, -t79 * t85; 0, t77, t77, 0, 0; 0, t76, t76, 0, 0; 0, t79, t79, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end