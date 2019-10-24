% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PPRRP2
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
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:19
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5PPRRP2_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP2_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP2_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PPRRP2_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPRRP2_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP2_jacobiaD_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:19:41
	% EndTime: 2019-10-24 10:19:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:19:41
	% EndTime: 2019-10-24 10:19:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:19:41
	% EndTime: 2019-10-24 10:19:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:19:41
	% EndTime: 2019-10-24 10:19:41
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->3), mult. (16->8), div. (0->0), fcn. (10->4), ass. (0->5)
	t14 = pkin(8) + qJ(3);
	t12 = sin(t14);
	t13 = cos(t14);
	t17 = qJD(3) * (-r_i_i_C(1) * t13 + r_i_i_C(2) * t12);
	t1 = [0, 0, cos(pkin(7)) * t17, 0, 0; 0, 0, sin(pkin(7)) * t17, 0, 0; 0, 0, (-r_i_i_C(1) * t12 - r_i_i_C(2) * t13) * qJD(3), 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:19:41
	% EndTime: 2019-10-24 10:19:41
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (64->20), mult. (116->41), div. (0->0), fcn. (87->6), ass. (0->20)
	t162 = pkin(6) + r_i_i_C(3);
	t147 = sin(pkin(7));
	t149 = sin(qJ(4));
	t161 = t147 * t149;
	t150 = cos(qJ(4));
	t160 = t147 * t150;
	t148 = cos(pkin(7));
	t159 = t148 * t149;
	t158 = t148 * t150;
	t146 = pkin(8) + qJ(3);
	t145 = cos(t146);
	t157 = qJD(3) * t145;
	t144 = sin(t146);
	t156 = qJD(4) * t144;
	t155 = r_i_i_C(1) * t149 + r_i_i_C(2) * t150;
	t154 = -r_i_i_C(1) * t150 + r_i_i_C(2) * t149 - pkin(3);
	t153 = t155 * t144;
	t152 = qJD(3) * t153;
	t151 = qJD(4) * t153 + (-t162 * t144 + t154 * t145) * qJD(3);
	t1 = [0, 0, t151 * t148, t148 * t152 + ((-t145 * t158 - t161) * r_i_i_C(1) + (t145 * t159 - t160) * r_i_i_C(2)) * qJD(4), 0; 0, 0, t151 * t147, t147 * t152 + ((-t145 * t160 + t159) * r_i_i_C(1) + (t145 * t161 + t158) * r_i_i_C(2)) * qJD(4), 0; 0, 0, -t155 * t145 * qJD(4) + (t154 * t144 + t162 * t145) * qJD(3), (t149 * t156 - t150 * t157) * r_i_i_C(2) + (-t149 * t157 - t150 * t156) * r_i_i_C(1), 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:19:41
	% EndTime: 2019-10-24 10:19:42
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (141->29), mult. (246->49), div. (0->0), fcn. (199->6), ass. (0->28)
	t200 = sin(qJ(4));
	t201 = cos(qJ(4));
	t217 = r_i_i_C(3) + qJ(5);
	t219 = pkin(4) + r_i_i_C(1);
	t220 = t219 * t200 - t217 * t201;
	t221 = t220 * qJD(4) - qJD(5) * t200;
	t218 = pkin(6) + r_i_i_C(2);
	t197 = pkin(8) + qJ(3);
	t195 = sin(t197);
	t216 = t195 * t201;
	t198 = sin(pkin(7));
	t215 = t198 * t200;
	t214 = t198 * t201;
	t199 = cos(pkin(7));
	t213 = t199 * t200;
	t212 = t199 * t201;
	t211 = qJD(3) * t200;
	t209 = t195 * t211;
	t208 = qJD(3) * t216;
	t196 = cos(t197);
	t207 = -t196 * t212 - t215;
	t206 = -t196 * t214 + t213;
	t205 = -t217 * t200 - t219 * t201;
	t203 = -pkin(3) + t205;
	t202 = t221 * t195 + (-t218 * t195 + t203 * t196) * qJD(3);
	t193 = t207 * qJD(4) + t199 * t209;
	t191 = t206 * qJD(4) + t198 * t209;
	t1 = [0, 0, t202 * t199, -t207 * qJD(5) - t217 * (t199 * t208 + (t196 * t213 - t214) * qJD(4)) + t219 * t193, -t193; 0, 0, t202 * t198, -t206 * qJD(5) - t217 * (t198 * t208 + (t196 * t215 + t212) * qJD(4)) + t219 * t191, -t191; 0, 0, -t221 * t196 + (t203 * t195 + t218 * t196) * qJD(3), -t220 * t196 * qJD(3) + (t205 * qJD(4) + qJD(5) * t201) * t195, qJD(4) * t216 + t196 * t211;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end