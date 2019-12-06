% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PRPRP1
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5PRPRP1_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP1_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP1_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PRPRP1_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPRP1_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP1_jacobiaD_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:29:18
	% EndTime: 2019-12-05 15:29:18
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:29:18
	% EndTime: 2019-12-05 15:29:18
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:29:18
	% EndTime: 2019-12-05 15:29:18
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (6->3), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->4)
	t32 = pkin(7) + qJ(2);
	t31 = cos(t32);
	t30 = sin(t32);
	t1 = [0, (-r_i_i_C(1) * t31 + r_i_i_C(2) * t30) * qJD(2), 0, 0, 0; 0, (-r_i_i_C(1) * t30 - r_i_i_C(2) * t31) * qJD(2), 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:29:18
	% EndTime: 2019-12-05 15:29:18
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (24->8), mult. (28->12), div. (0->0), fcn. (18->4), ass. (0->6)
	t19 = r_i_i_C(3) + qJ(3);
	t18 = -r_i_i_C(1) * cos(pkin(8)) + r_i_i_C(2) * sin(pkin(8)) - pkin(2);
	t15 = pkin(7) + qJ(2);
	t14 = cos(t15);
	t13 = sin(t15);
	t1 = [0, t14 * qJD(3) + (-t19 * t13 + t18 * t14) * qJD(2), qJD(2) * t14, 0, 0; 0, t13 * qJD(3) + (t18 * t13 + t19 * t14) * qJD(2), qJD(2) * t13, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:29:18
	% EndTime: 2019-12-05 15:29:18
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (69->20), mult. (70->32), div. (0->0), fcn. (46->5), ass. (0->15)
	t40 = r_i_i_C(3) + pkin(6) + qJ(3);
	t31 = pkin(7) + qJ(2);
	t27 = sin(t31);
	t39 = qJD(2) * t27;
	t29 = cos(t31);
	t38 = qJD(2) * t29;
	t37 = qJD(4) * t27;
	t36 = qJD(4) * t29;
	t30 = pkin(8) + qJ(4);
	t26 = sin(t30);
	t28 = cos(t30);
	t35 = r_i_i_C(1) * t26 + r_i_i_C(2) * t28;
	t34 = -r_i_i_C(1) * t28 + r_i_i_C(2) * t26 - cos(pkin(8)) * pkin(3) - pkin(2);
	t33 = t35 * qJD(4);
	t1 = [0, t29 * qJD(3) + t35 * t37 + (-t40 * t27 + t34 * t29) * qJD(2), t38, (t26 * t36 + t28 * t39) * r_i_i_C(2) + (t26 * t39 - t28 * t36) * r_i_i_C(1), 0; 0, t27 * qJD(3) - t29 * t33 + (t34 * t27 + t40 * t29) * qJD(2), t39, (t26 * t37 - t28 * t38) * r_i_i_C(2) + (-t26 * t38 - t28 * t37) * r_i_i_C(1), 0; 0, 0, 0, -t33, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:29:19
	% EndTime: 2019-12-05 15:29:19
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (146->25), mult. (140->37), div. (0->0), fcn. (98->5), ass. (0->17)
	t151 = pkin(8) + qJ(4);
	t147 = sin(t151);
	t149 = cos(t151);
	t162 = r_i_i_C(3) + qJ(5);
	t164 = pkin(4) + r_i_i_C(1);
	t165 = t164 * t147 - t162 * t149;
	t166 = t165 * qJD(4) - t147 * qJD(5);
	t163 = r_i_i_C(2) + pkin(6) + qJ(3);
	t152 = pkin(7) + qJ(2);
	t148 = sin(t152);
	t161 = qJD(2) * t148;
	t150 = cos(t152);
	t160 = qJD(2) * t150;
	t159 = qJD(4) * t150;
	t157 = -t162 * t147 - t164 * t149;
	t155 = -cos(pkin(8)) * pkin(3) - pkin(2) + t157;
	t1 = [0, t150 * qJD(3) + t166 * t148 + (-t163 * t148 + t155 * t150) * qJD(2), t160, (-t162 * t159 + t164 * t161) * t147 + (-t162 * t161 + (-t164 * qJD(4) + qJD(5)) * t150) * t149, -t147 * t161 + t149 * t159; 0, t148 * qJD(3) - t166 * t150 + (t155 * t148 + t163 * t150) * qJD(2), t161, -t165 * t160 + (t157 * qJD(4) + qJD(5) * t149) * t148, t148 * qJD(4) * t149 + t147 * t160; 0, 0, 0, -t166, qJD(4) * t147;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end