% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RPRRRP6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:54
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RPRRRP6_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP6_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP6_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRP6_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP6_jacobigD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:54:08
	% EndTime: 2019-10-10 01:54:08
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:54:08
	% EndTime: 2019-10-10 01:54:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:54:08
	% EndTime: 2019-10-10 01:54:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:54:08
	% EndTime: 2019-10-10 01:54:08
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (2->2), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, 0, qJD(1) * cos(qJ(1)), 0, 0, 0; 0, 0, qJD(1) * sin(qJ(1)), 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:54:08
	% EndTime: 2019-10-10 01:54:08
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (7->3), mult. (11->8), div. (0->0), fcn. (11->4), ass. (0->8)
	t82 = sin(qJ(1));
	t86 = qJD(1) * t82;
	t83 = cos(qJ(1));
	t85 = qJD(1) * t83;
	t81 = pkin(10) + qJ(3);
	t84 = qJD(3) * cos(t81);
	t79 = sin(t81);
	t1 = [0, 0, t85, -t79 * t86 + t83 * t84, 0, 0; 0, 0, t86, t79 * t85 + t82 * t84, 0, 0; 0, 0, 0, qJD(3) * t79, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:54:08
	% EndTime: 2019-10-10 01:54:08
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (14->3), mult. (20->8), div. (0->0), fcn. (20->4), ass. (0->11)
	t108 = sin(qJ(1));
	t112 = qJD(1) * t108;
	t109 = cos(qJ(1));
	t111 = qJD(1) * t109;
	t107 = pkin(10) + qJ(3);
	t110 = qJD(3) * cos(t107);
	t105 = sin(t107);
	t104 = qJD(3) * t105;
	t103 = t105 * t111 + t108 * t110;
	t102 = -t105 * t112 + t109 * t110;
	t1 = [0, 0, t111, t102, t102, 0; 0, 0, t112, t103, t103, 0; 0, 0, 0, t104, t104, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:54:08
	% EndTime: 2019-10-10 01:54:08
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (14->3), mult. (20->8), div. (0->0), fcn. (20->4), ass. (0->11)
	t109 = sin(qJ(1));
	t113 = qJD(1) * t109;
	t110 = cos(qJ(1));
	t112 = qJD(1) * t110;
	t108 = pkin(10) + qJ(3);
	t111 = qJD(3) * cos(t108);
	t106 = sin(t108);
	t105 = qJD(3) * t106;
	t104 = t106 * t112 + t109 * t111;
	t103 = -t106 * t113 + t110 * t111;
	t1 = [0, 0, t112, t103, t103, 0; 0, 0, t113, t104, t104, 0; 0, 0, 0, t105, t105, 0;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end