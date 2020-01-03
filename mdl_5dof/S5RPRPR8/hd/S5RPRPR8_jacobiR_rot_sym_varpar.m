% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPRPR8
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RPRPR8_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR8_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPR8_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR8_jacobiR_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:22:51
	% EndTime: 2019-12-31 18:22:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:22:51
	% EndTime: 2019-12-31 18:22:51
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
	% StartTime: 2019-12-31 18:22:51
	% EndTime: 2019-12-31 18:22:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (7->4), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->4)
	t12 = qJ(1) + pkin(8);
	t11 = cos(t12);
	t10 = sin(t12);
	t1 = [-t10, 0, 0, 0, 0; t11, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t11, 0, 0, 0, 0; -t10, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:22:51
	% EndTime: 2019-12-31 18:22:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (17->8), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->10)
	t13 = qJ(1) + pkin(8);
	t11 = sin(t13);
	t14 = sin(qJ(3));
	t19 = t11 * t14;
	t15 = cos(qJ(3));
	t18 = t11 * t15;
	t12 = cos(t13);
	t17 = t12 * t14;
	t16 = t12 * t15;
	t1 = [-t18, 0, -t17, 0, 0; t16, 0, -t19, 0, 0; 0, 0, t15, 0, 0; t19, 0, -t16, 0, 0; -t17, 0, -t18, 0, 0; 0, 0, -t14, 0, 0; t12, 0, 0, 0, 0; t11, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:22:51
	% EndTime: 2019-12-31 18:22:51
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (24->9), mult. (26->20), div. (0->0), fcn. (45->6), ass. (0->12)
	t52 = sin(pkin(9));
	t54 = sin(qJ(3));
	t59 = t52 * t54;
	t53 = cos(pkin(9));
	t58 = t53 * t54;
	t55 = cos(qJ(3));
	t57 = t55 * t52;
	t56 = t55 * t53;
	t51 = qJ(1) + pkin(8);
	t50 = cos(t51);
	t49 = sin(t51);
	t1 = [-t49 * t56 + t50 * t52, 0, -t50 * t58, 0, 0; t49 * t52 + t50 * t56, 0, -t49 * t58, 0, 0; 0, 0, t56, 0, 0; t49 * t57 + t50 * t53, 0, t50 * t59, 0, 0; t49 * t53 - t50 * t57, 0, t49 * t59, 0, 0; 0, 0, -t57, 0, 0; -t49 * t54, 0, t50 * t55, 0, 0; t50 * t54, 0, t49 * t55, 0, 0; 0, 0, t54, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:22:51
	% EndTime: 2019-12-31 18:22:51
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (62->14), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->17)
	t72 = pkin(9) + qJ(5);
	t68 = sin(t72);
	t74 = sin(qJ(3));
	t79 = t74 * t68;
	t70 = cos(t72);
	t78 = t74 * t70;
	t75 = cos(qJ(3));
	t77 = t75 * t68;
	t76 = t75 * t70;
	t73 = qJ(1) + pkin(8);
	t71 = cos(t73);
	t69 = sin(t73);
	t67 = t69 * t68 + t71 * t76;
	t66 = t69 * t70 - t71 * t77;
	t65 = t71 * t68 - t69 * t76;
	t64 = t69 * t77 + t71 * t70;
	t1 = [t65, 0, -t71 * t78, 0, t66; t67, 0, -t69 * t78, 0, -t64; 0, 0, t76, 0, -t79; t64, 0, t71 * t79, 0, -t67; t66, 0, t69 * t79, 0, t65; 0, 0, -t77, 0, -t78; -t69 * t74, 0, t71 * t75, 0, 0; t71 * t74, 0, t69 * t75, 0, 0; 0, 0, t74, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end