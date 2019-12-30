% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S5RRPRP10
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% JgD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 18:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S5RRPRP10_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP10_jacobigD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP10_jacobigD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRP10_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP10_jacobigD_rot_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:59:18
	% EndTime: 2019-12-29 18:59:18
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:59:18
	% EndTime: 2019-12-29 18:59:18
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:59:13
	% EndTime: 2019-12-29 18:59:13
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (2->2), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, qJD(1) * cos(qJ(1)), 0, 0, 0; 0, qJD(1) * sin(qJ(1)), 0, 0, 0; 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:59:18
	% EndTime: 2019-12-29 18:59:18
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (2->2), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, qJD(1) * cos(qJ(1)), 0, 0, 0; 0, qJD(1) * sin(qJ(1)), 0, 0, 0; 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:59:13
	% EndTime: 2019-12-29 18:59:13
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (2->2), mult. (11->8), div. (0->0), fcn. (11->4), ass. (0->7)
	t78 = sin(qJ(1));
	t83 = qJD(1) * t78;
	t80 = cos(qJ(1));
	t82 = qJD(1) * t80;
	t81 = qJD(2) * sin(qJ(2));
	t79 = cos(qJ(2));
	t1 = [0, t82, 0, -t79 * t83 - t80 * t81, 0; 0, t83, 0, -t78 * t81 + t79 * t82, 0; 0, 0, 0, qJD(2) * t79, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:59:13
	% EndTime: 2019-12-29 18:59:13
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (2->2), mult. (11->8), div. (0->0), fcn. (11->4), ass. (0->7)
	t82 = sin(qJ(1));
	t87 = qJD(1) * t82;
	t84 = cos(qJ(1));
	t86 = qJD(1) * t84;
	t85 = qJD(2) * sin(qJ(2));
	t83 = cos(qJ(2));
	t1 = [0, t86, 0, -t83 * t87 - t84 * t85, 0; 0, t87, 0, -t82 * t85 + t83 * t86, 0; 0, 0, 0, qJD(2) * t83, 0;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,5);
end