% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PRRPRP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:16
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRPRP1_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP1_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP1_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRP1_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP1_jacobigD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:16:24
	% EndTime: 2019-10-09 22:16:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:16:24
	% EndTime: 2019-10-09 22:16:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:16:24
	% EndTime: 2019-10-09 22:16:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:16:24
	% EndTime: 2019-10-09 22:16:24
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (2->2), mult. (12->9), div. (0->0), fcn. (12->6), ass. (0->6)
	t82 = sin(qJ(2));
	t84 = cos(pkin(6)) * t82;
	t83 = cos(qJ(2));
	t80 = cos(pkin(10));
	t79 = sin(pkin(10));
	t1 = [0, 0, (-t79 * t84 + t80 * t83) * qJD(2), 0, 0, 0; 0, 0, (t79 * t83 + t80 * t84) * qJD(2), 0, 0, 0; 0, 0, sin(pkin(6)) * qJD(2) * t82, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:16:24
	% EndTime: 2019-10-09 22:16:24
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (12->9), div. (0->0), fcn. (12->6), ass. (0->6)
	t90 = sin(qJ(2));
	t92 = cos(pkin(6)) * t90;
	t91 = cos(qJ(2));
	t88 = cos(pkin(10));
	t87 = sin(pkin(10));
	t1 = [0, 0, (-t87 * t92 + t88 * t91) * qJD(2), 0, 0, 0; 0, 0, (t87 * t91 + t88 * t92) * qJD(2), 0, 0, 0; 0, 0, sin(pkin(6)) * qJD(2) * t90, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:16:25
	% EndTime: 2019-10-09 22:16:25
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (21->11), mult. (48->30), div. (0->0), fcn. (50->8), ass. (0->17)
	t139 = qJ(3) + pkin(11);
	t137 = sin(t139);
	t141 = sin(pkin(6));
	t152 = t141 * t137;
	t143 = cos(pkin(6));
	t144 = sin(qJ(2));
	t151 = t143 * t144;
	t145 = cos(qJ(2));
	t150 = t143 * t145;
	t149 = qJD(2) * t137;
	t148 = qJD(2) * t141;
	t140 = sin(pkin(10));
	t142 = cos(pkin(10));
	t147 = t140 * t145 + t142 * t151;
	t146 = -t140 * t151 + t142 * t145;
	t138 = cos(t139);
	t1 = [0, 0, t146 * qJD(2), 0, (t146 * t138 + t140 * t152) * qJD(3) + (-t140 * t150 - t142 * t144) * t149, 0; 0, 0, t147 * qJD(2), 0, (t147 * t138 - t142 * t152) * qJD(3) + (-t140 * t144 + t142 * t150) * t149, 0; 0, 0, t144 * t148, 0, t145 * t137 * t148 + (t138 * t141 * t144 + t137 * t143) * qJD(3), 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:16:25
	% EndTime: 2019-10-09 22:16:25
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (21->11), mult. (48->30), div. (0->0), fcn. (50->8), ass. (0->17)
	t146 = qJ(3) + pkin(11);
	t144 = sin(t146);
	t148 = sin(pkin(6));
	t159 = t148 * t144;
	t150 = cos(pkin(6));
	t151 = sin(qJ(2));
	t158 = t150 * t151;
	t152 = cos(qJ(2));
	t157 = t150 * t152;
	t156 = qJD(2) * t144;
	t155 = qJD(2) * t148;
	t147 = sin(pkin(10));
	t149 = cos(pkin(10));
	t154 = t147 * t152 + t149 * t158;
	t153 = -t147 * t158 + t149 * t152;
	t145 = cos(t146);
	t1 = [0, 0, t153 * qJD(2), 0, (t153 * t145 + t147 * t159) * qJD(3) + (-t147 * t157 - t149 * t151) * t156, 0; 0, 0, t154 * qJD(2), 0, (t154 * t145 - t149 * t159) * qJD(3) + (-t147 * t151 + t149 * t157) * t156, 0; 0, 0, t151 * t155, 0, t152 * t144 * t155 + (t145 * t148 * t151 + t144 * t150) * qJD(3), 0;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end