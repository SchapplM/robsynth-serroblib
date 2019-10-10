% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PRPRPR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:30
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRPRPR1_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRPR1_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_jacobig_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:29:59
	% EndTime: 2019-10-09 21:29:59
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:29:59
	% EndTime: 2019-10-09 21:29:59
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:29:59
	% EndTime: 2019-10-09 21:29:59
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t18 = sin(pkin(6));
	t1 = [0, sin(pkin(10)) * t18, 0, 0, 0, 0; 0, -cos(pkin(10)) * t18, 0, 0, 0, 0; 0, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:29:59
	% EndTime: 2019-10-09 21:29:59
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t30 = sin(pkin(6));
	t1 = [0, sin(pkin(10)) * t30, 0, 0, 0, 0; 0, -cos(pkin(10)) * t30, 0, 0, 0, 0; 0, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:29:59
	% EndTime: 2019-10-09 21:29:59
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (8->5), mult. (22->12), div. (0->0), fcn. (35->8), ass. (0->12)
	t71 = sin(pkin(11));
	t74 = cos(pkin(11));
	t77 = sin(qJ(2));
	t78 = cos(qJ(2));
	t79 = t71 * t77 - t74 * t78;
	t76 = cos(pkin(6));
	t75 = cos(pkin(10));
	t73 = sin(pkin(6));
	t72 = sin(pkin(10));
	t70 = -t78 * t71 - t77 * t74;
	t69 = t79 * t76;
	t1 = [0, t72 * t73, 0, -t72 * t69 - t75 * t70, 0, 0; 0, -t75 * t73, 0, t75 * t69 - t72 * t70, 0, 0; 0, t76, 0, t79 * t73, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:29:59
	% EndTime: 2019-10-09 21:29:59
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (8->5), mult. (22->12), div. (0->0), fcn. (35->8), ass. (0->12)
	t80 = sin(pkin(11));
	t83 = cos(pkin(11));
	t86 = sin(qJ(2));
	t87 = cos(qJ(2));
	t88 = t80 * t86 - t83 * t87;
	t85 = cos(pkin(6));
	t84 = cos(pkin(10));
	t82 = sin(pkin(6));
	t81 = sin(pkin(10));
	t79 = -t87 * t80 - t86 * t83;
	t78 = t88 * t85;
	t1 = [0, t81 * t82, 0, -t81 * t78 - t84 * t79, 0, 0; 0, -t84 * t82, 0, t84 * t78 - t81 * t79, 0, 0; 0, t85, 0, t88 * t82, 0, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:29:59
	% EndTime: 2019-10-09 21:29:59
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (24->11), mult. (50->24), div. (0->0), fcn. (76->10), ass. (0->18)
	t122 = sin(pkin(10));
	t123 = sin(pkin(6));
	t132 = t122 * t123;
	t125 = cos(pkin(10));
	t131 = t125 * t123;
	t121 = sin(pkin(11));
	t124 = cos(pkin(11));
	t127 = sin(qJ(2));
	t128 = cos(qJ(2));
	t130 = t128 * t121 + t127 * t124;
	t129 = t127 * t121 - t128 * t124;
	t126 = cos(pkin(6));
	t120 = qJ(4) + pkin(12);
	t119 = cos(t120);
	t118 = sin(t120);
	t115 = t130 * t126;
	t114 = t129 * t126;
	t1 = [0, t132, 0, -t122 * t114 + t125 * t130, 0, (-t122 * t115 - t125 * t129) * t118 - t119 * t132; 0, -t131, 0, t125 * t114 + t122 * t130, 0, (t125 * t115 - t122 * t129) * t118 + t119 * t131; 0, t126, 0, t129 * t123, 0, t130 * t118 * t123 - t126 * t119;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end