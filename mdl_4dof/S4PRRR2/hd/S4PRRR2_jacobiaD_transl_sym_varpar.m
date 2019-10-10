% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S4PRRR2
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
% 
% Output:
% JaD_transl [3x4]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 20:37
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S4PRRR2_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(3,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR2_jacobiaD_transl_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR2_jacobiaD_transl_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S4PRRR2_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4PRRR2_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S4PRRR2_jacobiaD_transl_sym_varpar: pkin has to be [2x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:37:38
	% EndTime: 2019-10-09 20:37:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:37:38
	% EndTime: 2019-10-09 20:37:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:37:38
	% EndTime: 2019-10-09 20:37:38
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t6 = cos(qJ(2));
	t5 = sin(qJ(2));
	t1 = [0, (r_i_i_C(1) * t5 + r_i_i_C(2) * t6) * qJD(2), 0, 0; 0, 0, 0, 0; 0, (-r_i_i_C(1) * t6 + r_i_i_C(2) * t5) * qJD(2), 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:37:38
	% EndTime: 2019-10-09 20:37:38
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (22->6), mult. (20->9), div. (0->0), fcn. (10->4), ass. (0->8)
	t15 = qJD(2) + qJD(3);
	t16 = qJ(2) + qJ(3);
	t21 = sin(t16) * t15;
	t20 = cos(t16) * t15;
	t19 = r_i_i_C(1) * t21 + r_i_i_C(2) * t20;
	t18 = pkin(1) * qJD(2);
	t17 = -r_i_i_C(1) * t20 + r_i_i_C(2) * t21;
	t1 = [0, sin(qJ(2)) * t18 + t19, t19, 0; 0, 0, 0, 0; 0, -cos(qJ(2)) * t18 + t17, t17, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:37:38
	% EndTime: 2019-10-09 20:37:38
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (68->10), mult. (36->12), div. (0->0), fcn. (18->6), ass. (0->13)
	t24 = qJD(2) + qJD(3);
	t33 = pkin(2) * t24;
	t21 = qJD(4) + t24;
	t25 = qJ(2) + qJ(3);
	t23 = qJ(4) + t25;
	t32 = sin(t23) * t21;
	t31 = cos(t23) * t21;
	t30 = r_i_i_C(1) * t32 + r_i_i_C(2) * t31;
	t29 = pkin(1) * qJD(2);
	t28 = sin(t25) * t33 + t30;
	t27 = -r_i_i_C(1) * t31 + r_i_i_C(2) * t32;
	t26 = -cos(t25) * t33 + t27;
	t1 = [0, sin(qJ(2)) * t29 + t28, t28, t30; 0, 0, 0, 0; 0, -cos(qJ(2)) * t29 + t26, t26, t27;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,4);
end