% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRPPR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4,theta5]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:24
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRPPR6_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR6_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPPR6_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR6_jacobiR_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:24:00
	% EndTime: 2019-10-10 00:24:00
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:24:00
	% EndTime: 2019-10-10 00:24:00
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
	% StartTime: 2019-10-10 00:24:00
	% EndTime: 2019-10-10 00:24:00
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t6 = cos(qJ(1));
	t5 = sin(qJ(1));
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t5, 0, 0, 0, 0, 0; -t6, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t6, 0, 0, 0, 0, 0; t5, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:24:00
	% EndTime: 2019-10-10 00:24:00
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (5->5), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->9)
	t11 = sin(qJ(3));
	t12 = sin(qJ(1));
	t16 = t12 * t11;
	t13 = cos(qJ(3));
	t14 = cos(qJ(1));
	t15 = t14 * t13;
	t10 = t14 * t11;
	t9 = t12 * t13;
	t1 = [t10, 0, t9, 0, 0, 0; t16, 0, -t15, 0, 0, 0; 0, 0, -t11, 0, 0, 0; t15, 0, -t16, 0, 0, 0; t9, 0, t10, 0, 0, 0; 0, 0, -t13, 0, 0, 0; -t12, 0, 0, 0, 0, 0; t14, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:24:00
	% EndTime: 2019-10-10 00:24:00
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (15->6), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->10)
	t17 = qJ(3) + pkin(9);
	t15 = sin(t17);
	t18 = sin(qJ(1));
	t21 = t18 * t15;
	t16 = cos(t17);
	t19 = cos(qJ(1));
	t20 = t19 * t16;
	t14 = t19 * t15;
	t13 = t18 * t16;
	t1 = [t14, 0, t13, 0, 0, 0; t21, 0, -t20, 0, 0, 0; 0, 0, -t15, 0, 0, 0; t20, 0, -t21, 0, 0, 0; t13, 0, t14, 0, 0, 0; 0, 0, -t16, 0, 0, 0; -t18, 0, 0, 0, 0, 0; t19, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:24:01
	% EndTime: 2019-10-10 00:24:01
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (25->11), mult. (26->18), div. (0->0), fcn. (45->6), ass. (0->12)
	t57 = sin(pkin(10));
	t59 = sin(qJ(1));
	t64 = t59 * t57;
	t58 = cos(pkin(10));
	t63 = t59 * t58;
	t60 = cos(qJ(1));
	t62 = t60 * t57;
	t61 = t60 * t58;
	t56 = qJ(3) + pkin(9);
	t55 = cos(t56);
	t54 = sin(t56);
	t1 = [t54 * t61 - t64, 0, t55 * t63, 0, 0, 0; t54 * t63 + t62, 0, -t55 * t61, 0, 0, 0; 0, 0, -t54 * t58, 0, 0, 0; -t54 * t62 - t63, 0, -t55 * t64, 0, 0, 0; -t54 * t64 + t61, 0, t55 * t62, 0, 0, 0; 0, 0, t54 * t57, 0, 0, 0; -t60 * t55, 0, t59 * t54, 0, 0, 0; -t59 * t55, 0, -t60 * t54, 0, 0, 0; 0, 0, t55, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:24:01
	% EndTime: 2019-10-10 00:24:01
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (61->16), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->19)
	t81 = qJ(3) + pkin(9);
	t77 = sin(t81);
	t82 = sin(qJ(1));
	t89 = t82 * t77;
	t80 = pkin(10) + qJ(6);
	t78 = cos(t80);
	t88 = t82 * t78;
	t79 = cos(t81);
	t87 = t82 * t79;
	t83 = cos(qJ(1));
	t86 = t83 * t77;
	t85 = t83 * t78;
	t84 = t83 * t79;
	t76 = sin(t80);
	t75 = -t82 * t76 + t77 * t85;
	t74 = t76 * t86 + t88;
	t73 = t83 * t76 + t77 * t88;
	t72 = -t76 * t89 + t85;
	t1 = [t75, 0, t78 * t87, 0, 0, t72; t73, 0, -t78 * t84, 0, 0, t74; 0, 0, -t77 * t78, 0, 0, -t79 * t76; -t74, 0, -t76 * t87, 0, 0, -t73; t72, 0, t76 * t84, 0, 0, t75; 0, 0, t77 * t76, 0, 0, -t79 * t78; -t84, 0, t89, 0, 0, 0; -t87, 0, -t86, 0, 0, 0; 0, 0, t79, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end