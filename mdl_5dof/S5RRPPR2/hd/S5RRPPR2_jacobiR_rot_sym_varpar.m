% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRPPR2
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
%   Siehe auch: S5RRPPR2_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:06
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RRPPR2_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPPR2_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_jacobiR_rot_sym_varpar: pkin has to be [9x1] (double)');
JR_rot=NaN(9,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 10:06:21
	% EndTime: 2022-01-20 10:06:21
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 10:06:21
	% EndTime: 2022-01-20 10:06:21
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t8, 0, 0, 0, 0; t9, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t9, 0, 0, 0, 0; -t8, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 10:06:21
	% EndTime: 2022-01-20 10:06:21
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
	% StartTime: 2022-01-20 10:06:21
	% EndTime: 2022-01-20 10:06:21
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (22->7), mult. (0->0), div. (0->0), fcn. (8->2), ass. (0->4)
	t19 = qJ(1) + qJ(2) + pkin(8);
	t18 = cos(t19);
	t17 = sin(t19);
	t1 = [-t17, -t17, 0, 0, 0; t18, t18, 0, 0, 0; 0, 0, 0, 0, 0; -t18, -t18, 0, 0, 0; -t17, -t17, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 10:06:21
	% EndTime: 2022-01-20 10:06:22
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (28->5), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->10)
	t26 = qJ(1) + qJ(2) + pkin(8);
	t24 = sin(t26);
	t28 = cos(pkin(9));
	t30 = t24 * t28;
	t25 = cos(t26);
	t27 = sin(pkin(9));
	t29 = t25 * t27;
	t23 = t25 * t28;
	t22 = t24 * t27;
	t1 = [-t30, -t30, 0, 0, 0; t23, t23, 0, 0, 0; 0, 0, 0, 0, 0; t22, t22, 0, 0, 0; -t29, -t29, 0, 0, 0; 0, 0, 0, 0, 0; t25, t25, 0, 0, 0; t24, t24, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 10:06:22
	% EndTime: 2022-01-20 10:06:22
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (72->11), mult. (42->14), div. (0->0), fcn. (72->6), ass. (0->16)
	t72 = qJ(1) + qJ(2) + pkin(8);
	t70 = sin(t72);
	t73 = sin(pkin(9));
	t79 = t70 * t73;
	t74 = cos(pkin(9));
	t75 = sin(qJ(5));
	t78 = t74 * t75;
	t76 = cos(qJ(5));
	t77 = t74 * t76;
	t71 = cos(t72);
	t69 = t71 * t73;
	t68 = t70 * t75 + t71 * t77;
	t67 = t70 * t76 - t71 * t78;
	t66 = -t70 * t77 + t71 * t75;
	t65 = t70 * t78 + t71 * t76;
	t1 = [t66, t66, 0, 0, t67; t68, t68, 0, 0, -t65; 0, 0, 0, 0, -t73 * t75; t65, t65, 0, 0, -t68; t67, t67, 0, 0, t66; 0, 0, 0, 0, -t73 * t76; -t79, -t79, 0, 0, 0; t69, t69, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
end