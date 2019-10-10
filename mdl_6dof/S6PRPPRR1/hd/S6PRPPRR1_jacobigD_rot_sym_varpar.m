% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PRPPRR1
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3,theta4]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:24
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRPPRR1_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR1_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR1_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPPRR1_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPPRR1_jacobigD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:24:27
	% EndTime: 2019-10-09 21:24:27
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:24:27
	% EndTime: 2019-10-09 21:24:27
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:24:27
	% EndTime: 2019-10-09 21:24:27
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:24:27
	% EndTime: 2019-10-09 21:24:27
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:24:28
	% EndTime: 2019-10-09 21:24:28
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:24:28
	% EndTime: 2019-10-09 21:24:28
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (7->4), mult. (30->12), div. (0->0), fcn. (30->8), ass. (0->10)
	t114 = sin(pkin(11));
	t117 = cos(pkin(11));
	t120 = sin(qJ(2));
	t121 = cos(qJ(2));
	t123 = qJD(2) * (t114 * t121 + t117 * t120);
	t118 = cos(pkin(10));
	t115 = sin(pkin(10));
	t113 = (t114 * t120 - t117 * t121) * qJD(2);
	t112 = cos(pkin(6)) * t123;
	t1 = [0, 0, 0, 0, -t115 * t112 - t118 * t113, 0; 0, 0, 0, 0, t118 * t112 - t115 * t113, 0; 0, 0, 0, 0, sin(pkin(6)) * t123, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:24:28
	% EndTime: 2019-10-09 21:24:28
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (36->15), mult. (97->37), div. (0->0), fcn. (104->10), ass. (0->20)
	t182 = sin(pkin(11));
	t185 = cos(pkin(11));
	t188 = sin(qJ(2));
	t189 = cos(qJ(2));
	t191 = t188 * t182 - t189 * t185;
	t194 = t191 * qJD(2);
	t192 = t182 * t189 + t185 * t188;
	t177 = t192 * qJD(2);
	t181 = pkin(12) + qJ(5);
	t179 = sin(t181);
	t184 = sin(pkin(6));
	t193 = t184 * t179;
	t187 = cos(pkin(6));
	t186 = cos(pkin(10));
	t183 = sin(pkin(10));
	t180 = cos(t181);
	t175 = t192 * t187;
	t174 = t187 * t194;
	t173 = t187 * t177;
	t1 = [0, 0, 0, 0, -t183 * t173 - t186 * t194, (t183 * t174 - t186 * t177) * t179 + ((-t183 * t175 - t186 * t191) * t180 + t183 * t193) * qJD(5); 0, 0, 0, 0, t186 * t173 - t183 * t194, (-t186 * t174 - t183 * t177) * t179 + ((t186 * t175 - t183 * t191) * t180 - t186 * t193) * qJD(5); 0, 0, 0, 0, t184 * t177, t187 * qJD(5) * t179 + (t192 * qJD(5) * t180 - t179 * t194) * t184;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end