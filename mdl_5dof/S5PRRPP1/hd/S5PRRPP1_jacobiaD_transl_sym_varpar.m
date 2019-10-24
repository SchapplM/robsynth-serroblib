% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PRRPP1
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
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:28
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5PRRPP1_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP1_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP1_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PRRPP1_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRPP1_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP1_jacobiaD_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:28:45
	% EndTime: 2019-10-24 10:28:45
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:28:45
	% EndTime: 2019-10-24 10:28:45
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:28:45
	% EndTime: 2019-10-24 10:28:45
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
	% StartTime: 2019-10-24 10:28:45
	% EndTime: 2019-10-24 10:28:45
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (41->16), mult. (64->29), div. (0->0), fcn. (40->4), ass. (0->13)
	t24 = sin(qJ(3));
	t25 = cos(qJ(3));
	t26 = (r_i_i_C(1) * t24 + r_i_i_C(2) * t25) * qJD(3);
	t33 = pkin(6) + r_i_i_C(3);
	t32 = qJD(2) * t24;
	t31 = qJD(2) * t25;
	t30 = qJD(3) * t24;
	t29 = qJD(3) * t25;
	t27 = -r_i_i_C(1) * t25 + r_i_i_C(2) * t24 - pkin(2);
	t23 = pkin(7) + qJ(2);
	t22 = cos(t23);
	t21 = sin(t23);
	t1 = [0, t21 * t26 + (-t21 * t33 + t22 * t27) * qJD(2), (t21 * t31 + t22 * t30) * r_i_i_C(2) + (t21 * t32 - t22 * t29) * r_i_i_C(1), 0, 0; 0, -t22 * t26 + (t21 * t27 + t22 * t33) * qJD(2), (t21 * t30 - t22 * t31) * r_i_i_C(2) + (-t21 * t29 - t22 * t32) * r_i_i_C(1), 0, 0; 0, 0, -t26, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:28:45
	% EndTime: 2019-10-24 10:28:46
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (82->17), mult. (90->25), div. (0->0), fcn. (59->6), ass. (0->15)
	t34 = qJ(3) + pkin(8);
	t30 = sin(t34);
	t32 = cos(t34);
	t47 = -r_i_i_C(1) * t32 + r_i_i_C(2) * t30 - cos(qJ(3)) * pkin(3);
	t45 = r_i_i_C(3) + qJ(4) + pkin(6);
	t33 = pkin(7) + qJ(2);
	t31 = cos(t33);
	t44 = qJD(2) * t31;
	t42 = -pkin(2) + t47;
	t41 = sin(qJ(3)) * pkin(3) + r_i_i_C(1) * t30 + r_i_i_C(2) * t32;
	t29 = sin(t33);
	t40 = t41 * t29;
	t39 = qJD(3) * t47;
	t38 = t41 * qJD(3);
	t1 = [0, t31 * qJD(4) + qJD(3) * t40 + (-t45 * t29 + t42 * t31) * qJD(2), qJD(2) * t40 + t31 * t39, t44, 0; 0, t29 * qJD(4) - t31 * t38 + (t42 * t29 + t45 * t31) * qJD(2), t29 * t39 - t41 * t44, qJD(2) * t29, 0; 0, 0, -t38, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:28:46
	% EndTime: 2019-10-24 10:28:46
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (159->21), mult. (160->28), div. (0->0), fcn. (111->6), ass. (0->18)
	t155 = qJ(3) + pkin(8);
	t151 = sin(t155);
	t153 = cos(t155);
	t169 = r_i_i_C(3) + qJ(5);
	t173 = pkin(4) + r_i_i_C(1);
	t162 = t173 * t151 - t169 * t153 + sin(qJ(3)) * pkin(3);
	t159 = -t162 * qJD(3) + t151 * qJD(5);
	t176 = t159 + (r_i_i_C(2) + qJ(4) + pkin(6)) * qJD(2);
	t175 = -t169 * t151 - t173 * t153 - cos(qJ(3)) * pkin(3);
	t154 = pkin(7) + qJ(2);
	t150 = sin(t154);
	t168 = qJD(2) * t150;
	t152 = cos(t154);
	t167 = qJD(2) * t152;
	t166 = qJD(3) * t153;
	t161 = qJD(4) + (-pkin(2) + t175) * qJD(2);
	t160 = t175 * qJD(3) + qJD(5) * t153;
	t1 = [0, -t176 * t150 + t161 * t152, t160 * t152 + t162 * t168, t167, -t151 * t168 + t152 * t166; 0, t161 * t150 + t176 * t152, t160 * t150 - t162 * t167, t168, t150 * t166 + t151 * t167; 0, 0, t159, 0, qJD(3) * t151;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end