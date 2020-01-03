% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPPPR6
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RPPPR6_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR6_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPPR6_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR6_jacobiR_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:48:16
	% EndTime: 2019-12-31 17:48:16
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:48:16
	% EndTime: 2019-12-31 17:48:16
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
	% StartTime: 2019-12-31 17:48:16
	% EndTime: 2019-12-31 17:48:16
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (10->4), ass. (0->5)
	t10 = cos(qJ(1));
	t9 = sin(qJ(1));
	t8 = cos(pkin(7));
	t7 = sin(pkin(7));
	t1 = [-t9 * t8, 0, 0, 0, 0; t10 * t8, 0, 0, 0, 0; 0, 0, 0, 0, 0; t9 * t7, 0, 0, 0, 0; -t10 * t7, 0, 0, 0, 0; 0, 0, 0, 0, 0; t10, 0, 0, 0, 0; t9, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:48:16
	% EndTime: 2019-12-31 17:48:16
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (10->4), ass. (0->5)
	t34 = cos(qJ(1));
	t33 = sin(qJ(1));
	t32 = cos(pkin(7));
	t31 = sin(pkin(7));
	t1 = [t34, 0, 0, 0, 0; t33, 0, 0, 0, 0; 0, 0, 0, 0, 0; t33 * t32, 0, 0, 0, 0; -t34 * t32, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t33 * t31, 0, 0, 0, 0; t34 * t31, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:48:16
	% EndTime: 2019-12-31 17:48:16
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (5->5), mult. (14->10), div. (0->0), fcn. (24->6), ass. (0->11)
	t37 = sin(pkin(8));
	t41 = sin(qJ(1));
	t46 = t41 * t37;
	t39 = cos(pkin(8));
	t45 = t41 * t39;
	t42 = cos(qJ(1));
	t44 = t42 * t37;
	t43 = t42 * t39;
	t40 = cos(pkin(7));
	t38 = sin(pkin(7));
	t1 = [-t38 * t46 + t43, 0, 0, 0, 0; t38 * t44 + t45, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t38 * t45 - t44, 0, 0, 0, 0; t38 * t43 - t46, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t41 * t40, 0, 0, 0, 0; t42 * t40, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:48:16
	% EndTime: 2019-12-31 17:48:16
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (20->12), mult. (60->23), div. (0->0), fcn. (92->8), ass. (0->22)
	t77 = cos(pkin(7));
	t78 = sin(qJ(5));
	t90 = t77 * t78;
	t80 = cos(qJ(5));
	t89 = t77 * t80;
	t81 = cos(qJ(1));
	t88 = t77 * t81;
	t74 = sin(pkin(8));
	t79 = sin(qJ(1));
	t87 = t79 * t74;
	t76 = cos(pkin(8));
	t86 = t79 * t76;
	t85 = t81 * t74;
	t84 = t81 * t76;
	t75 = sin(pkin(7));
	t72 = -t75 * t87 + t84;
	t83 = t72 * t78 + t79 * t89;
	t82 = t72 * t80 - t79 * t90;
	t71 = t75 * t85 + t86;
	t70 = t71 * t80 + t78 * t88;
	t69 = -t71 * t78 + t80 * t88;
	t1 = [t82, 0, 0, 0, t69; t70, 0, 0, 0, t83; 0, 0, 0, 0, t74 * t90 + t75 * t80; -t83, 0, 0, 0, -t70; t69, 0, 0, 0, t82; 0, 0, 0, 0, t74 * t89 - t75 * t78; t75 * t86 + t85, 0, 0, 0, 0; -t75 * t84 + t87, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end