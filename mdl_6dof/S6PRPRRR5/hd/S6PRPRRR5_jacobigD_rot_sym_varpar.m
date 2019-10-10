% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PRPRRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:01
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRPRRR5_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR5_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR5_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRRR5_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR5_jacobigD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:01:16
	% EndTime: 2019-10-09 22:01:16
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:01:16
	% EndTime: 2019-10-09 22:01:16
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:01:16
	% EndTime: 2019-10-09 22:01:16
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:01:16
	% EndTime: 2019-10-09 22:01:16
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:01:16
	% EndTime: 2019-10-09 22:01:17
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (12->9), div. (0->0), fcn. (12->6), ass. (0->6)
	t82 = cos(qJ(2));
	t83 = cos(pkin(6)) * t82;
	t81 = sin(qJ(2));
	t79 = cos(pkin(11));
	t78 = sin(pkin(11));
	t1 = [0, 0, 0, (-t78 * t83 - t79 * t81) * qJD(2), 0, 0; 0, 0, 0, (-t78 * t81 + t79 * t83) * qJD(2), 0, 0; 0, 0, 0, sin(pkin(6)) * qJD(2) * t82, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:01:17
	% EndTime: 2019-10-09 22:01:17
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (4->2), mult. (24->9), div. (0->0), fcn. (24->6), ass. (0->9)
	t111 = cos(qJ(2));
	t112 = cos(pkin(6)) * t111;
	t110 = sin(qJ(2));
	t108 = cos(pkin(11));
	t107 = sin(pkin(11));
	t106 = sin(pkin(6)) * qJD(2) * t111;
	t105 = (-t107 * t112 - t108 * t110) * qJD(2);
	t104 = (-t107 * t110 + t108 * t112) * qJD(2);
	t1 = [0, 0, 0, t105, t105, 0; 0, 0, 0, t104, t104, 0; 0, 0, 0, t106, t106, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:01:17
	% EndTime: 2019-10-09 22:01:17
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (29->14), mult. (60->30), div. (0->0), fcn. (62->8), ass. (0->21)
	t166 = qJ(4) + qJ(5);
	t164 = cos(t166);
	t168 = sin(pkin(6));
	t179 = t168 * t164;
	t170 = cos(pkin(6));
	t171 = sin(qJ(2));
	t178 = t170 * t171;
	t172 = cos(qJ(2));
	t177 = t170 * t172;
	t176 = qJD(2) * t164;
	t175 = qJD(2) * t168;
	t167 = sin(pkin(11));
	t169 = cos(pkin(11));
	t174 = -t167 * t171 + t169 * t177;
	t173 = t167 * t177 + t169 * t171;
	t165 = qJD(4) + qJD(5);
	t163 = sin(t166);
	t162 = t172 * t175;
	t161 = t173 * qJD(2);
	t160 = t174 * qJD(2);
	t1 = [0, 0, 0, -t161, -t161, (t173 * t163 + t167 * t179) * t165 - (-t167 * t178 + t169 * t172) * t176; 0, 0, 0, t160, t160, (-t174 * t163 - t169 * t179) * t165 - (t167 * t172 + t169 * t178) * t176; 0, 0, 0, t162, t162, -t168 * t172 * t165 * t163 + (t165 * t170 - t171 * t175) * t164;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end