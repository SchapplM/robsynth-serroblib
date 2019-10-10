% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S4RPPP1
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
% 
% Output:
% JRD_rot [9x4]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 20:39
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S4RPPP1_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_jacobiRD_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPP1_jacobiRD_rot_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4RPPP1_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_jacobiRD_rot_sym_varpar: pkin has to be [6x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:39:01
	% EndTime: 2019-10-09 20:39:01
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:39:01
	% EndTime: 2019-10-09 20:39:01
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0; -t31, 0, 0, 0; 0, 0, 0, 0; t31, 0, 0, 0; -t30, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:39:02
	% EndTime: 2019-10-09 20:39:02
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (5->5), mult. (24->17), div. (0->0), fcn. (24->6), ass. (0->9)
	t128 = cos(pkin(4));
	t129 = sin(qJ(1));
	t133 = t128 * t129;
	t130 = cos(qJ(1));
	t132 = t128 * t130;
	t131 = qJD(1) * sin(pkin(4));
	t127 = cos(pkin(6));
	t125 = sin(pkin(6));
	t1 = [(t125 * t133 - t127 * t130) * qJD(1), 0, 0, 0; (-t125 * t132 - t127 * t129) * qJD(1), 0, 0, 0; 0, 0, 0, 0; (t125 * t130 + t127 * t133) * qJD(1), 0, 0, 0; (t125 * t129 - t127 * t132) * qJD(1), 0, 0, 0; 0, 0, 0, 0; -t129 * t131, 0, 0, 0; t130 * t131, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:39:02
	% EndTime: 2019-10-09 20:39:02
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (5->5), mult. (24->17), div. (0->0), fcn. (24->6), ass. (0->9)
	t162 = cos(pkin(4));
	t163 = sin(qJ(1));
	t167 = t162 * t163;
	t164 = cos(qJ(1));
	t166 = t162 * t164;
	t165 = qJD(1) * sin(pkin(4));
	t161 = cos(pkin(6));
	t159 = sin(pkin(6));
	t1 = [-t163 * t165, 0, 0, 0; t164 * t165, 0, 0, 0; 0, 0, 0, 0; (-t159 * t167 + t161 * t164) * qJD(1), 0, 0, 0; (t159 * t166 + t161 * t163) * qJD(1), 0, 0, 0; 0, 0, 0, 0; (-t159 * t164 - t161 * t167) * qJD(1), 0, 0, 0; (-t159 * t163 + t161 * t166) * qJD(1), 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:39:02
	% EndTime: 2019-10-09 20:39:02
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (5->5), mult. (24->17), div. (0->0), fcn. (24->6), ass. (0->9)
	t164 = cos(pkin(4));
	t165 = sin(qJ(1));
	t169 = t164 * t165;
	t166 = cos(qJ(1));
	t168 = t164 * t166;
	t167 = qJD(1) * sin(pkin(4));
	t163 = cos(pkin(6));
	t161 = sin(pkin(6));
	t1 = [-t165 * t167, 0, 0, 0; t166 * t167, 0, 0, 0; 0, 0, 0, 0; (-t161 * t166 - t163 * t169) * qJD(1), 0, 0, 0; (-t161 * t165 + t163 * t168) * qJD(1), 0, 0, 0; 0, 0, 0, 0; (t161 * t169 - t163 * t166) * qJD(1), 0, 0, 0; (-t161 * t168 - t163 * t165) * qJD(1), 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,4);
end