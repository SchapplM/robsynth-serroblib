% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PRRRPR2
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:48
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRRRPR2_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR2_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRPR2_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR2_jacobig_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:48:33
	% EndTime: 2019-10-09 22:48:33
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:48:33
	% EndTime: 2019-10-09 22:48:33
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:48:33
	% EndTime: 2019-10-09 22:48:33
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t18 = sin(pkin(6));
	t1 = [0, sin(pkin(11)) * t18, 0, 0, 0, 0; 0, -cos(pkin(11)) * t18, 0, 0, 0, 0; 0, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:48:34
	% EndTime: 2019-10-09 22:48:34
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (4->4), mult. (9->8), div. (0->0), fcn. (17->6), ass. (0->8)
	t57 = cos(pkin(6));
	t59 = cos(qJ(2));
	t60 = t57 * t59;
	t58 = sin(qJ(2));
	t56 = cos(pkin(11));
	t55 = sin(pkin(6));
	t54 = sin(pkin(11));
	t1 = [0, t54 * t55, t54 * t60 + t56 * t58, 0, 0, 0; 0, -t56 * t55, t54 * t58 - t56 * t60, 0, 0, 0; 0, t57, -t55 * t59, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:48:34
	% EndTime: 2019-10-09 22:48:34
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (7->5), mult. (16->8), div. (0->0), fcn. (29->6), ass. (0->11)
	t73 = sin(pkin(6));
	t77 = cos(qJ(2));
	t79 = t73 * t77;
	t75 = cos(pkin(6));
	t78 = t75 * t77;
	t76 = sin(qJ(2));
	t74 = cos(pkin(11));
	t72 = sin(pkin(11));
	t71 = t72 * t78 + t74 * t76;
	t70 = t72 * t76 - t74 * t78;
	t1 = [0, t72 * t73, t71, t71, 0, 0; 0, -t74 * t73, t70, t70, 0, 0; 0, t75, -t79, -t79, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:48:34
	% EndTime: 2019-10-09 22:48:34
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (7->5), mult. (16->8), div. (0->0), fcn. (29->6), ass. (0->11)
	t110 = sin(pkin(6));
	t114 = cos(qJ(2));
	t116 = t110 * t114;
	t112 = cos(pkin(6));
	t115 = t112 * t114;
	t113 = sin(qJ(2));
	t111 = cos(pkin(11));
	t109 = sin(pkin(11));
	t108 = t109 * t115 + t111 * t113;
	t107 = t109 * t113 - t111 * t115;
	t1 = [0, t109 * t110, t108, t108, 0, 0; 0, -t111 * t110, t107, t107, 0, 0; 0, t112, -t116, -t116, 0, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:48:34
	% EndTime: 2019-10-09 22:48:34
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (18->11), mult. (31->20), div. (0->0), fcn. (52->8), ass. (0->17)
	t124 = sin(pkin(11));
	t125 = sin(pkin(6));
	t134 = t124 * t125;
	t129 = cos(qJ(2));
	t133 = t125 * t129;
	t126 = cos(pkin(11));
	t132 = t126 * t125;
	t127 = cos(pkin(6));
	t128 = sin(qJ(2));
	t131 = t127 * t128;
	t130 = t127 * t129;
	t123 = qJ(3) + qJ(4);
	t122 = cos(t123);
	t121 = sin(t123);
	t120 = t124 * t130 + t126 * t128;
	t119 = t124 * t128 - t126 * t130;
	t1 = [0, t134, t120, t120, 0, (-t124 * t131 + t126 * t129) * t121 - t122 * t134; 0, -t132, t119, t119, 0, (t124 * t129 + t126 * t131) * t121 + t122 * t132; 0, t127, -t133, -t133, 0, t125 * t128 * t121 - t127 * t122;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end