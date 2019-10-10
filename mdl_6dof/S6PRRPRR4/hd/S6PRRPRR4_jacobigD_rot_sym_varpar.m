% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PRRPRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:31
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRPRR4_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR4_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRR4_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_jacobigD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:31:29
	% EndTime: 2019-10-09 22:31:29
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:31:29
	% EndTime: 2019-10-09 22:31:29
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:31:29
	% EndTime: 2019-10-09 22:31:29
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:31:29
	% EndTime: 2019-10-09 22:31:29
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (12->9), div. (0->0), fcn. (12->6), ass. (0->6)
	t82 = sin(qJ(2));
	t84 = cos(pkin(6)) * t82;
	t83 = cos(qJ(2));
	t80 = cos(pkin(11));
	t79 = sin(pkin(11));
	t1 = [0, 0, (-t79 * t84 + t80 * t83) * qJD(2), 0, 0, 0; 0, 0, (t79 * t83 + t80 * t84) * qJD(2), 0, 0, 0; 0, 0, sin(pkin(6)) * qJD(2) * t82, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:31:29
	% EndTime: 2019-10-09 22:31:29
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (2->2), mult. (12->9), div. (0->0), fcn. (12->6), ass. (0->6)
	t100 = sin(qJ(2));
	t102 = t100 * cos(pkin(6));
	t101 = cos(qJ(2));
	t98 = cos(pkin(11));
	t97 = sin(pkin(11));
	t1 = [0, 0, (t101 * t98 - t97 * t102) * qJD(2), 0, 0, 0; 0, 0, (t101 * t97 + t98 * t102) * qJD(2), 0, 0, 0; 0, 0, sin(pkin(6)) * qJD(2) * t100, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:31:29
	% EndTime: 2019-10-09 22:31:29
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (5->5), mult. (24->9), div. (0->0), fcn. (24->6), ass. (0->9)
	t124 = sin(qJ(2));
	t127 = cos(pkin(6)) * t124;
	t126 = sin(pkin(6)) * qJD(2) * t124;
	t125 = cos(qJ(2));
	t122 = cos(pkin(11));
	t120 = sin(pkin(11));
	t119 = (-t120 * t127 + t122 * t125) * qJD(2);
	t118 = (t120 * t125 + t122 * t127) * qJD(2);
	t1 = [0, 0, t119, 0, -t119, 0; 0, 0, t118, 0, -t118, 0; 0, 0, t126, 0, -t126, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:31:30
	% EndTime: 2019-10-09 22:31:30
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (44->21), mult. (144->44), div. (0->0), fcn. (158->10), ass. (0->23)
	t250 = qJD(5) - qJD(3);
	t225 = sin(pkin(6));
	t229 = sin(qJ(3));
	t246 = t225 * t229;
	t230 = sin(qJ(2));
	t245 = t225 * t230;
	t232 = cos(qJ(3));
	t244 = t225 * t232;
	t227 = cos(pkin(6));
	t243 = t227 * t230;
	t233 = cos(qJ(2));
	t242 = t227 * t233;
	t241 = qJD(2) * t245;
	t224 = sin(pkin(11));
	t226 = cos(pkin(11));
	t222 = t224 * t233 + t226 * t243;
	t223 = -t224 * t243 + t226 * t233;
	t228 = sin(qJ(5));
	t231 = cos(qJ(5));
	t237 = qJD(2) * (t232 * t228 - t229 * t231);
	t221 = t223 * qJD(2);
	t219 = t222 * qJD(2);
	t1 = [0, 0, t221, 0, -t221, (-t224 * t242 - t226 * t230) * t237 - t250 * ((-t223 * t229 + t224 * t244) * t228 - (t223 * t232 + t224 * t246) * t231); 0, 0, t219, 0, -t219, (-t224 * t230 + t226 * t242) * t237 + t250 * ((t222 * t229 + t226 * t244) * t228 + (t222 * t232 - t226 * t246) * t231); 0, 0, t241, 0, -t241, t233 * t225 * t237 + t250 * ((-t227 * t232 + t229 * t245) * t228 + (t227 * t229 + t230 * t244) * t231);];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end