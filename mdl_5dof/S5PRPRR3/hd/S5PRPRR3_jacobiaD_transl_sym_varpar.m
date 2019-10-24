% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PRPRR3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:26
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5PRPRR3_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR3_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR3_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PRPRR3_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPRR3_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR3_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:26:20
	% EndTime: 2019-10-24 10:26:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:26:20
	% EndTime: 2019-10-24 10:26:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:26:20
	% EndTime: 2019-10-24 10:26:20
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->2), mult. (16->8), div. (0->0), fcn. (10->4), ass. (0->4)
	t12 = sin(qJ(2));
	t13 = cos(qJ(2));
	t14 = qJD(2) * (-r_i_i_C(1) * t13 + r_i_i_C(2) * t12);
	t1 = [0, cos(pkin(8)) * t14, 0, 0, 0; 0, sin(pkin(8)) * t14, 0, 0, 0; 0, (-r_i_i_C(1) * t12 - r_i_i_C(2) * t13) * qJD(2), 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:26:20
	% EndTime: 2019-10-24 10:26:20
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (12->5), mult. (24->10), div. (0->0), fcn. (15->6), ass. (0->5)
	t15 = qJ(2) + pkin(9);
	t13 = sin(t15);
	t14 = cos(t15);
	t19 = qJD(2) * (-cos(qJ(2)) * pkin(2) - r_i_i_C(1) * t14 + r_i_i_C(2) * t13);
	t1 = [0, cos(pkin(8)) * t19, 0, 0, 0; 0, sin(pkin(8)) * t19, 0, 0, 0; 0, (-sin(qJ(2)) * pkin(2) - r_i_i_C(1) * t13 - r_i_i_C(2) * t14) * qJD(2), 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:26:21
	% EndTime: 2019-10-24 10:26:21
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (67->22), mult. (124->43), div. (0->0), fcn. (92->8), ass. (0->20)
	t164 = pkin(6) + r_i_i_C(3);
	t148 = sin(pkin(8));
	t150 = sin(qJ(4));
	t163 = t148 * t150;
	t151 = cos(qJ(4));
	t162 = t148 * t151;
	t149 = cos(pkin(8));
	t161 = t149 * t150;
	t160 = t149 * t151;
	t147 = qJ(2) + pkin(9);
	t146 = cos(t147);
	t159 = qJD(2) * t146;
	t145 = sin(t147);
	t158 = qJD(4) * t145;
	t157 = r_i_i_C(1) * t150 + r_i_i_C(2) * t151;
	t156 = -r_i_i_C(1) * t151 + r_i_i_C(2) * t150 - pkin(3);
	t155 = t157 * t145;
	t154 = qJD(2) * t155;
	t153 = qJD(4) * t155 + (-cos(qJ(2)) * pkin(2) - t164 * t145 + t156 * t146) * qJD(2);
	t1 = [0, t153 * t149, 0, t149 * t154 + ((-t146 * t160 - t163) * r_i_i_C(1) + (t146 * t161 - t162) * r_i_i_C(2)) * qJD(4), 0; 0, t153 * t148, 0, t148 * t154 + ((-t146 * t162 + t161) * r_i_i_C(1) + (t146 * t163 + t160) * r_i_i_C(2)) * qJD(4), 0; 0, -t157 * t146 * qJD(4) + (-sin(qJ(2)) * pkin(2) + t164 * t146 + t156 * t145) * qJD(2), 0, (t150 * t158 - t151 * t159) * r_i_i_C(2) + (-t150 * t159 - t151 * t158) * r_i_i_C(1), 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:26:21
	% EndTime: 2019-10-24 10:26:21
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (188->36), mult. (214->63), div. (0->0), fcn. (163->10), ass. (0->34)
	t186 = qJD(4) + qJD(5);
	t191 = sin(qJ(4));
	t188 = qJ(4) + qJ(5);
	t184 = sin(t188);
	t185 = cos(t188);
	t214 = r_i_i_C(2) * t185;
	t199 = r_i_i_C(1) * t184 + t214;
	t212 = pkin(4) * qJD(4);
	t216 = t186 * t199 + t191 * t212;
	t215 = r_i_i_C(2) * t184;
	t213 = r_i_i_C(3) + pkin(7) + pkin(6);
	t187 = qJ(2) + pkin(9);
	t183 = cos(t187);
	t192 = cos(qJ(4));
	t211 = t183 * t192;
	t210 = t185 * t186;
	t189 = sin(pkin(8));
	t209 = t186 * t189;
	t190 = cos(pkin(8));
	t208 = t186 * t190;
	t182 = sin(t187);
	t205 = qJD(2) * t182;
	t196 = t189 * t205 + t208;
	t202 = t183 * t209;
	t207 = (t184 * t196 - t185 * t202) * r_i_i_C(1) + (t184 * t202 + t185 * t196) * r_i_i_C(2);
	t197 = t190 * t205 - t209;
	t201 = t183 * t208;
	t206 = (t184 * t197 - t185 * t201) * r_i_i_C(1) + (t184 * t201 + t185 * t197) * r_i_i_C(2);
	t204 = qJD(2) * t183;
	t200 = t191 * t205;
	t198 = -t192 * pkin(4) - r_i_i_C(1) * t185 - pkin(3) + t215;
	t195 = t216 * t182 + (-cos(qJ(2)) * pkin(2) - t213 * t182 + t198 * t183) * qJD(2);
	t180 = t182 * t186 * t215;
	t1 = [0, t195 * t190, 0, (t190 * t200 + (-t189 * t191 - t190 * t211) * qJD(4)) * pkin(4) + t206, t206; 0, t195 * t189, 0, (t189 * t200 + (-t189 * t211 + t190 * t191) * qJD(4)) * pkin(4) + t207, t207; 0, -t216 * t183 + (-sin(qJ(2)) * pkin(2) + t213 * t183 + t198 * t182) * qJD(2), 0, t180 + (-r_i_i_C(1) * t210 - t192 * t212) * t182 + (-pkin(4) * t191 - t199) * t204, -t204 * t214 + t180 + (-t182 * t210 - t184 * t204) * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end