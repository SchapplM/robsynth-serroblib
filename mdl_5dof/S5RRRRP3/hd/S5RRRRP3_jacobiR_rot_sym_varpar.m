% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRRRP3
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RRRRP3_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP3_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRP3_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP3_jacobiR_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:49:50
	% EndTime: 2019-12-31 21:49:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:49:51
	% EndTime: 2019-12-31 21:49:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t8, 0, 0, 0, 0; t9, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t9, 0, 0, 0, 0; -t8, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:49:50
	% EndTime: 2019-12-31 21:49:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (14->7), mult. (0->0), div. (0->0), fcn. (8->2), ass. (0->4)
	t16 = qJ(1) + qJ(2);
	t15 = cos(t16);
	t14 = sin(t16);
	t1 = [-t14, -t14, 0, 0, 0; t15, t15, 0, 0, 0; 0, 0, 0, 0, 0; -t15, -t15, 0, 0, 0; -t14, -t14, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:49:51
	% EndTime: 2019-12-31 21:49:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (33->10), mult. (0->0), div. (0->0), fcn. (12->2), ass. (0->4)
	t21 = qJ(1) + qJ(2) + qJ(3);
	t20 = cos(t21);
	t19 = sin(t21);
	t1 = [-t19, -t19, -t19, 0, 0; t20, t20, t20, 0, 0; 0, 0, 0, 0, 0; -t20, -t20, -t20, 0, 0; -t19, -t19, -t19, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:49:51
	% EndTime: 2019-12-31 21:49:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (55->12), mult. (16->4), div. (0->0), fcn. (40->4), ass. (0->10)
	t31 = qJ(1) + qJ(2) + qJ(3);
	t29 = sin(t31);
	t33 = cos(qJ(4));
	t35 = t29 * t33;
	t30 = cos(t31);
	t32 = sin(qJ(4));
	t34 = t30 * t32;
	t28 = t30 * t33;
	t27 = t29 * t32;
	t1 = [-t35, -t35, -t35, -t34, 0; t28, t28, t28, -t27, 0; 0, 0, 0, t33, 0; t27, t27, t27, -t28, 0; -t34, -t34, -t34, -t35, 0; 0, 0, 0, -t32, 0; t30, t30, t30, 0, 0; t29, t29, t29, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:49:51
	% EndTime: 2019-12-31 21:49:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (52->9), mult. (16->4), div. (0->0), fcn. (40->4), ass. (0->10)
	t65 = qJ(1) + qJ(2) + qJ(3);
	t63 = sin(t65);
	t66 = sin(qJ(4));
	t69 = t63 * t66;
	t67 = cos(qJ(4));
	t68 = t63 * t67;
	t64 = cos(t65);
	t62 = t64 * t67;
	t61 = t64 * t66;
	t1 = [-t68, -t68, -t68, -t61, 0; t62, t62, t62, -t69, 0; 0, 0, 0, t67, 0; t64, t64, t64, 0, 0; t63, t63, t63, 0, 0; 0, 0, 0, 0, 0; -t69, -t69, -t69, t62, 0; t61, t61, t61, t68, 0; 0, 0, 0, t66, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end