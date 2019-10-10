% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RPRPRR8
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:58
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RPRPRR8_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR8_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR8_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR8_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR8_jacobigD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:58:22
	% EndTime: 2019-10-10 00:58:22
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:58:22
	% EndTime: 2019-10-10 00:58:22
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:58:22
	% EndTime: 2019-10-10 00:58:22
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:58:22
	% EndTime: 2019-10-10 00:58:22
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, 0, -qJD(1) * sin(qJ(1)), 0, 0, 0; 0, 0, qJD(1) * cos(qJ(1)), 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:58:22
	% EndTime: 2019-10-10 00:58:22
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, 0, -qJD(1) * sin(qJ(1)), 0, 0, 0; 0, 0, qJD(1) * cos(qJ(1)), 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:58:22
	% EndTime: 2019-10-10 00:58:22
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (8->4), mult. (11->8), div. (0->0), fcn. (11->4), ass. (0->8)
	t85 = sin(qJ(1));
	t89 = qJD(1) * t85;
	t86 = cos(qJ(1));
	t88 = qJD(1) * t86;
	t84 = qJ(3) + pkin(10);
	t87 = qJD(3) * sin(t84);
	t83 = cos(t84);
	t1 = [0, 0, -t89, 0, -t83 * t88 + t85 * t87, 0; 0, 0, t88, 0, -t83 * t89 - t86 * t87, 0; 0, 0, 0, 0, qJD(3) * t83, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:58:22
	% EndTime: 2019-10-10 00:58:22
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (15->4), mult. (20->8), div. (0->0), fcn. (20->4), ass. (0->11)
	t110 = sin(qJ(1));
	t114 = qJD(1) * t110;
	t111 = cos(qJ(1));
	t113 = qJD(1) * t111;
	t109 = qJ(3) + pkin(10);
	t112 = qJD(3) * sin(t109);
	t108 = cos(t109);
	t106 = qJD(3) * t108;
	t105 = -t108 * t113 + t110 * t112;
	t104 = -t108 * t114 - t111 * t112;
	t1 = [0, 0, -t114, 0, t105, t105; 0, 0, t113, 0, t104, t104; 0, 0, 0, 0, t106, t106;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end