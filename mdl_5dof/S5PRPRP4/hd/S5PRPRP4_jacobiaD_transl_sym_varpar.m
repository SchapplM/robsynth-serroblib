% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PRPRP4
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
% Datum: 2019-12-05 15:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5PRPRP4_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP4_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP4_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PRPRP4_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPRP4_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP4_jacobiaD_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:36:31
	% EndTime: 2019-12-05 15:36:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:36:31
	% EndTime: 2019-12-05 15:36:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:36:31
	% EndTime: 2019-12-05 15:36:31
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->2), mult. (16->8), div. (0->0), fcn. (10->4), ass. (0->4)
	t12 = sin(qJ(2));
	t13 = cos(qJ(2));
	t14 = qJD(2) * (-r_i_i_C(1) * t13 + r_i_i_C(2) * t12);
	t1 = [0, cos(pkin(7)) * t14, 0, 0, 0; 0, sin(pkin(7)) * t14, 0, 0, 0; 0, (-r_i_i_C(1) * t12 - r_i_i_C(2) * t13) * qJD(2), 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:36:31
	% EndTime: 2019-12-05 15:36:31
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (12->5), mult. (24->10), div. (0->0), fcn. (15->6), ass. (0->5)
	t15 = qJ(2) + pkin(8);
	t13 = sin(t15);
	t14 = cos(t15);
	t19 = qJD(2) * (-cos(qJ(2)) * pkin(2) - r_i_i_C(1) * t14 + r_i_i_C(2) * t13);
	t1 = [0, cos(pkin(7)) * t19, 0, 0, 0; 0, sin(pkin(7)) * t19, 0, 0, 0; 0, (-sin(qJ(2)) * pkin(2) - r_i_i_C(1) * t13 - r_i_i_C(2) * t14) * qJD(2), 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:36:31
	% EndTime: 2019-12-05 15:36:31
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (67->22), mult. (124->43), div. (0->0), fcn. (92->8), ass. (0->20)
	t164 = pkin(6) + r_i_i_C(3);
	t148 = sin(pkin(7));
	t150 = sin(qJ(4));
	t163 = t148 * t150;
	t151 = cos(qJ(4));
	t162 = t148 * t151;
	t149 = cos(pkin(7));
	t161 = t149 * t150;
	t160 = t149 * t151;
	t147 = qJ(2) + pkin(8);
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
	% StartTime: 2019-12-05 15:36:32
	% EndTime: 2019-12-05 15:36:32
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (144->31), mult. (254->51), div. (0->0), fcn. (204->8), ass. (0->28)
	t201 = sin(qJ(4));
	t202 = cos(qJ(4));
	t219 = r_i_i_C(3) + qJ(5);
	t221 = pkin(4) + r_i_i_C(1);
	t222 = t221 * t201 - t219 * t202;
	t223 = t222 * qJD(4) - qJD(5) * t201;
	t220 = pkin(6) + r_i_i_C(2);
	t198 = qJ(2) + pkin(8);
	t196 = sin(t198);
	t218 = t196 * t202;
	t199 = sin(pkin(7));
	t217 = t199 * t201;
	t216 = t199 * t202;
	t200 = cos(pkin(7));
	t215 = t200 * t201;
	t214 = t200 * t202;
	t213 = qJD(2) * t201;
	t211 = t196 * t213;
	t210 = qJD(2) * t218;
	t197 = cos(t198);
	t209 = -t197 * t214 - t217;
	t208 = -t197 * t216 + t215;
	t207 = -t219 * t201 - t221 * t202;
	t205 = -pkin(3) + t207;
	t204 = t223 * t196 + (-cos(qJ(2)) * pkin(2) - t220 * t196 + t205 * t197) * qJD(2);
	t194 = t209 * qJD(4) + t200 * t211;
	t192 = t208 * qJD(4) + t199 * t211;
	t1 = [0, t204 * t200, 0, -t209 * qJD(5) - t219 * (t200 * t210 + (t197 * t215 - t216) * qJD(4)) + t221 * t194, -t194; 0, t204 * t199, 0, -t208 * qJD(5) - t219 * (t199 * t210 + (t197 * t217 + t214) * qJD(4)) + t221 * t192, -t192; 0, -t223 * t197 + (-sin(qJ(2)) * pkin(2) + t220 * t197 + t205 * t196) * qJD(2), 0, -t222 * t197 * qJD(2) + (t207 * qJD(4) + qJD(5) * t202) * t196, qJD(4) * t218 + t197 * t213;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end