% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PRPRPR5
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:37
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRPRPR5_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR5_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR5_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRPR5_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR5_jacobigD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:37:21
	% EndTime: 2019-10-09 21:37:21
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:37:21
	% EndTime: 2019-10-09 21:37:21
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:37:21
	% EndTime: 2019-10-09 21:37:21
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:37:21
	% EndTime: 2019-10-09 21:37:21
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:37:21
	% EndTime: 2019-10-09 21:37:21
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (12->9), div. (0->0), fcn. (12->6), ass. (0->6)
	t87 = sin(qJ(2));
	t89 = cos(pkin(6)) * t87;
	t88 = cos(qJ(2));
	t85 = cos(pkin(10));
	t84 = sin(pkin(10));
	t1 = [0, 0, 0, (-t84 * t89 + t85 * t88) * qJD(2), 0, 0; 0, 0, 0, (t84 * t88 + t85 * t89) * qJD(2), 0, 0; 0, 0, 0, sin(pkin(6)) * qJD(2) * t87, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:37:21
	% EndTime: 2019-10-09 21:37:21
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (12->9), div. (0->0), fcn. (12->6), ass. (0->6)
	t106 = sin(qJ(2));
	t108 = cos(pkin(6)) * t106;
	t107 = cos(qJ(2));
	t104 = cos(pkin(10));
	t103 = sin(pkin(10));
	t1 = [0, 0, 0, (-t103 * t108 + t104 * t107) * qJD(2), 0, 0; 0, 0, 0, (t103 * t107 + t104 * t108) * qJD(2), 0, 0; 0, 0, 0, sin(pkin(6)) * qJD(2) * t106, 0, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:37:22
	% EndTime: 2019-10-09 21:37:22
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (21->11), mult. (48->30), div. (0->0), fcn. (50->8), ass. (0->17)
	t137 = pkin(11) + qJ(4);
	t136 = cos(t137);
	t139 = sin(pkin(6));
	t150 = t139 * t136;
	t141 = cos(pkin(6));
	t142 = sin(qJ(2));
	t149 = t141 * t142;
	t143 = cos(qJ(2));
	t148 = t141 * t143;
	t147 = qJD(2) * t136;
	t146 = qJD(2) * t139;
	t138 = sin(pkin(10));
	t140 = cos(pkin(10));
	t145 = t138 * t143 + t140 * t149;
	t144 = -t138 * t149 + t140 * t143;
	t135 = sin(t137);
	t1 = [0, 0, 0, t144 * qJD(2), 0, (-t144 * t135 + t138 * t150) * qJD(4) + (-t138 * t148 - t140 * t142) * t147; 0, 0, 0, t145 * qJD(2), 0, (-t145 * t135 - t140 * t150) * qJD(4) + (-t138 * t142 + t140 * t148) * t147; 0, 0, 0, t142 * t146, 0, t143 * t136 * t146 + (-t135 * t139 * t142 + t136 * t141) * qJD(4);];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end