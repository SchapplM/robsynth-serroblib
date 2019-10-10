% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PRPRRR6
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:03
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRPRRR6_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR6_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRRR6_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR6_jacobig_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:03:06
	% EndTime: 2019-10-09 22:03:07
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:03:06
	% EndTime: 2019-10-09 22:03:07
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:03:07
	% EndTime: 2019-10-09 22:03:07
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t18 = sin(pkin(6));
	t1 = [0, sin(pkin(11)) * t18, 0, 0, 0, 0; 0, -cos(pkin(11)) * t18, 0, 0, 0, 0; 0, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:03:07
	% EndTime: 2019-10-09 22:03:07
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t42 = sin(pkin(6));
	t1 = [0, sin(pkin(11)) * t42, 0, 0, 0, 0; 0, -cos(pkin(11)) * t42, 0, 0, 0, 0; 0, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:03:07
	% EndTime: 2019-10-09 22:03:07
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->3), mult. (9->8), div. (0->0), fcn. (17->6), ass. (0->8)
	t57 = cos(pkin(6));
	t58 = sin(qJ(2));
	t60 = t57 * t58;
	t59 = cos(qJ(2));
	t56 = cos(pkin(11));
	t55 = sin(pkin(6));
	t54 = sin(pkin(11));
	t1 = [0, t54 * t55, 0, -t54 * t60 + t56 * t59, 0, 0; 0, -t56 * t55, 0, t54 * t59 + t56 * t60, 0, 0; 0, t57, 0, t55 * t58, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:03:07
	% EndTime: 2019-10-09 22:03:07
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (8->8), mult. (24->20), div. (0->0), fcn. (40->8), ass. (0->13)
	t84 = sin(pkin(11));
	t85 = sin(pkin(6));
	t95 = t84 * t85;
	t86 = cos(pkin(11));
	t94 = t86 * t85;
	t87 = cos(pkin(6));
	t89 = sin(qJ(2));
	t93 = t87 * t89;
	t91 = cos(qJ(2));
	t92 = t87 * t91;
	t90 = cos(qJ(4));
	t88 = sin(qJ(4));
	t1 = [0, t95, 0, -t84 * t93 + t86 * t91, t88 * t95 - (t84 * t92 + t86 * t89) * t90, 0; 0, -t94, 0, t84 * t91 + t86 * t93, -t88 * t94 - (t84 * t89 - t86 * t92) * t90, 0; 0, t87, 0, t85 * t89, t85 * t90 * t91 + t87 * t88, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:03:07
	% EndTime: 2019-10-09 22:03:07
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (13->8), mult. (39->20), div. (0->0), fcn. (63->8), ass. (0->16)
	t102 = sin(pkin(11));
	t103 = sin(pkin(6));
	t113 = t102 * t103;
	t104 = cos(pkin(11));
	t112 = t104 * t103;
	t105 = cos(pkin(6));
	t107 = sin(qJ(2));
	t111 = t105 * t107;
	t109 = cos(qJ(2));
	t110 = t105 * t109;
	t108 = cos(qJ(4));
	t106 = sin(qJ(4));
	t101 = t103 * t109 * t108 + t105 * t106;
	t100 = -t106 * t112 - (t102 * t107 - t104 * t110) * t108;
	t99 = t106 * t113 - (t102 * t110 + t104 * t107) * t108;
	t1 = [0, t113, 0, -t102 * t111 + t104 * t109, t99, t99; 0, -t112, 0, t102 * t109 + t104 * t111, t100, t100; 0, t105, 0, t103 * t107, t101, t101;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end