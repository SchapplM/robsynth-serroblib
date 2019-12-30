% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRPPR6
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 18:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RRPPR6_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR6_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPPR6_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR6_jacobiR_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:19:54
	% EndTime: 2019-12-29 18:19:54
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:19:49
	% EndTime: 2019-12-29 18:19:49
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t8, 0, 0, 0, 0; t9, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t9, 0, 0, 0, 0; -t8, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:19:49
	% EndTime: 2019-12-29 18:19:49
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (7->7), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->9)
	t10 = sin(qJ(1));
	t9 = sin(qJ(2));
	t16 = t10 * t9;
	t12 = cos(qJ(1));
	t15 = t12 * t9;
	t11 = cos(qJ(2));
	t14 = t10 * t11;
	t13 = t12 * t11;
	t1 = [-t14, -t15, 0, 0, 0; t13, -t16, 0, 0, 0; 0, t11, 0, 0, 0; t16, -t13, 0, 0, 0; -t15, -t14, 0, 0, 0; 0, -t9, 0, 0, 0; t12, 0, 0, 0, 0; t10, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:19:49
	% EndTime: 2019-12-29 18:19:49
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (17->8), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->10)
	t17 = qJ(2) + pkin(8);
	t15 = sin(t17);
	t18 = sin(qJ(1));
	t23 = t18 * t15;
	t16 = cos(t17);
	t22 = t18 * t16;
	t19 = cos(qJ(1));
	t21 = t19 * t15;
	t20 = t19 * t16;
	t1 = [-t22, -t21, 0, 0, 0; t20, -t23, 0, 0, 0; 0, t16, 0, 0, 0; t23, -t20, 0, 0, 0; -t21, -t22, 0, 0, 0; 0, -t15, 0, 0, 0; t19, 0, 0, 0, 0; t18, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:19:54
	% EndTime: 2019-12-29 18:19:54
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (23->9), mult. (26->18), div. (0->0), fcn. (45->6), ass. (0->12)
	t57 = sin(pkin(9));
	t59 = sin(qJ(1));
	t64 = t59 * t57;
	t58 = cos(pkin(9));
	t63 = t59 * t58;
	t60 = cos(qJ(1));
	t62 = t60 * t57;
	t61 = t60 * t58;
	t56 = qJ(2) + pkin(8);
	t55 = cos(t56);
	t54 = sin(t56);
	t1 = [-t55 * t63 + t62, -t54 * t61, 0, 0, 0; t55 * t61 + t64, -t54 * t63, 0, 0, 0; 0, t55 * t58, 0, 0, 0; t55 * t64 + t61, t54 * t62, 0, 0, 0; -t55 * t62 + t63, t54 * t64, 0, 0, 0; 0, -t55 * t57, 0, 0, 0; -t59 * t54, t60 * t55, 0, 0, 0; t60 * t54, t59 * t55, 0, 0, 0; 0, t54, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:19:49
	% EndTime: 2019-12-29 18:19:49
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (59->14), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->19)
	t80 = qJ(2) + pkin(8);
	t76 = sin(t80);
	t81 = sin(qJ(1));
	t88 = t81 * t76;
	t79 = pkin(9) + qJ(5);
	t77 = cos(t79);
	t87 = t81 * t77;
	t78 = cos(t80);
	t86 = t81 * t78;
	t82 = cos(qJ(1));
	t85 = t82 * t76;
	t84 = t82 * t77;
	t83 = t82 * t78;
	t75 = sin(t79);
	t74 = t81 * t75 + t77 * t83;
	t73 = -t75 * t83 + t87;
	t72 = t82 * t75 - t77 * t86;
	t71 = t75 * t86 + t84;
	t1 = [t72, -t76 * t84, 0, 0, t73; t74, -t76 * t87, 0, 0, -t71; 0, t78 * t77, 0, 0, -t76 * t75; t71, t75 * t85, 0, 0, -t74; t73, t75 * t88, 0, 0, t72; 0, -t78 * t75, 0, 0, -t76 * t77; -t88, t83, 0, 0, 0; t85, t86, 0, 0, 0; 0, t76, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end