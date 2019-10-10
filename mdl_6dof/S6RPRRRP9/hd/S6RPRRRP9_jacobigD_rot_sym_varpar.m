% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RPRRRP9
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 08:52
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RPRRRP9_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP9_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP9_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRP9_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP9_jacobigD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:52:22
	% EndTime: 2019-10-10 08:52:22
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:52:22
	% EndTime: 2019-10-10 08:52:22
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:52:22
	% EndTime: 2019-10-10 08:52:22
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:52:22
	% EndTime: 2019-10-10 08:52:22
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, 0, -qJD(1) * sin(qJ(1)), 0, 0, 0; 0, 0, qJD(1) * cos(qJ(1)), 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:52:22
	% EndTime: 2019-10-10 08:52:22
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (3->3), mult. (11->8), div. (0->0), fcn. (11->4), ass. (0->7)
	t78 = sin(qJ(1));
	t83 = qJD(1) * t78;
	t80 = cos(qJ(1));
	t82 = qJD(1) * t80;
	t81 = qJD(3) * sin(qJ(3));
	t79 = cos(qJ(3));
	t1 = [0, 0, -t83, t78 * t81 - t79 * t82, 0, 0; 0, 0, t82, -t79 * t83 - t80 * t81, 0, 0; 0, 0, 0, qJD(3) * t79, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:52:22
	% EndTime: 2019-10-10 08:52:22
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (5->3), mult. (20->8), div. (0->0), fcn. (20->4), ass. (0->10)
	t96 = sin(qJ(1));
	t101 = qJD(1) * t96;
	t98 = cos(qJ(1));
	t100 = qJD(1) * t98;
	t99 = qJD(3) * sin(qJ(3));
	t97 = cos(qJ(3));
	t94 = qJD(3) * t97;
	t93 = -t97 * t100 + t96 * t99;
	t92 = -t97 * t101 - t98 * t99;
	t1 = [0, 0, -t101, t93, t93, 0; 0, 0, t100, t92, t92, 0; 0, 0, 0, t94, t94, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:52:22
	% EndTime: 2019-10-10 08:52:22
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (5->3), mult. (20->8), div. (0->0), fcn. (20->4), ass. (0->10)
	t98 = sin(qJ(1));
	t103 = qJD(1) * t98;
	t102 = qJD(3) * sin(qJ(3));
	t100 = cos(qJ(1));
	t101 = qJD(1) * t100;
	t99 = cos(qJ(3));
	t96 = qJD(3) * t99;
	t95 = -t99 * t101 + t98 * t102;
	t94 = -t100 * t102 - t99 * t103;
	t1 = [0, 0, -t103, t95, t95, 0; 0, 0, t101, t94, t94, 0; 0, 0, 0, t96, t96, 0;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end