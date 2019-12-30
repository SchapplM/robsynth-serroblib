% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S4RPPR4
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
% 
% Output:
% JaD_transl [3x4]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 12:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S4RPPR4_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR4_jacobiaD_transl_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR4_jacobiaD_transl_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S4RPPR4_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4RPPR4_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR4_jacobiaD_transl_sym_varpar: pkin has to be [6x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:40:12
	% EndTime: 2019-12-29 12:40:12
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:40:12
	% EndTime: 2019-12-29 12:40:12
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0; 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:40:12
	% EndTime: 2019-12-29 12:40:12
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (8->5), mult. (12->8), div. (0->0), fcn. (6->4), ass. (0->4)
	t32 = qJ(1) + pkin(6);
	t31 = cos(t32);
	t30 = sin(t32);
	t1 = [(-cos(qJ(1)) * pkin(1) - r_i_i_C(1) * t31 + r_i_i_C(2) * t30) * qJD(1), 0, 0, 0; (-sin(qJ(1)) * pkin(1) - r_i_i_C(1) * t30 - r_i_i_C(2) * t31) * qJD(1), 0, 0, 0; 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:40:12
	% EndTime: 2019-12-29 12:40:12
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (22->9), mult. (24->12), div. (0->0), fcn. (14->4), ass. (0->6)
	t13 = -pkin(2) + r_i_i_C(2);
	t12 = r_i_i_C(3) + qJ(3);
	t11 = qJ(1) + pkin(6);
	t10 = cos(t11);
	t9 = sin(t11);
	t1 = [t10 * qJD(3) + (-cos(qJ(1)) * pkin(1) - t12 * t9 + t13 * t10) * qJD(1), 0, qJD(1) * t10, 0; t9 * qJD(3) + (-sin(qJ(1)) * pkin(1) + t13 * t9 + t12 * t10) * qJD(1), 0, qJD(1) * t9, 0; 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:40:12
	% EndTime: 2019-12-29 12:40:13
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (53->19), mult. (76->33), div. (0->0), fcn. (48->6), ass. (0->14)
	t24 = sin(qJ(4));
	t25 = cos(qJ(4));
	t34 = (r_i_i_C(1) * t25 - r_i_i_C(2) * t24) * qJD(4);
	t33 = qJD(1) * t24;
	t32 = qJD(1) * t25;
	t31 = qJD(4) * t24;
	t30 = qJD(4) * t25;
	t29 = -pkin(2) - pkin(5) - r_i_i_C(3);
	t27 = r_i_i_C(1) * t24 + r_i_i_C(2) * t25 + qJ(3);
	t26 = qJD(3) + t34;
	t23 = qJ(1) + pkin(6);
	t22 = cos(t23);
	t21 = sin(t23);
	t1 = [t26 * t22 + (-cos(qJ(1)) * pkin(1) + t29 * t22 - t27 * t21) * qJD(1), 0, qJD(1) * t22, (-t21 * t30 - t22 * t33) * r_i_i_C(2) + (-t21 * t31 + t22 * t32) * r_i_i_C(1); t26 * t21 + (-sin(qJ(1)) * pkin(1) + t27 * t22 + t29 * t21) * qJD(1), 0, qJD(1) * t21, (-t21 * t33 + t22 * t30) * r_i_i_C(2) + (t21 * t32 + t22 * t31) * r_i_i_C(1); 0, 0, 0, -t34;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,4);
end