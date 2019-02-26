% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPR5
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:03
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRPR5_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR5_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR5_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPR5_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR5_jacobiaD_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:03:13
% EndTime: 2019-02-26 21:03:13
% DurationCPUTime: 0.18s
% Computational Cost: add. (263->42), mult. (223->50), div. (0->0), fcn. (154->7), ass. (0->36)
t193 = pkin(10) + qJ(3);
t191 = qJ(4) + t193;
t188 = cos(t191);
t220 = r_i_i_C(3) + qJ(5);
t204 = t220 * t188;
t187 = sin(t191);
t186 = t187 * qJD(5);
t194 = qJD(3) + qJD(4);
t189 = sin(t193);
t219 = pkin(3) * qJD(3);
t212 = t189 * t219;
t223 = pkin(4) - r_i_i_C(2);
t228 = (-t223 * t187 + t204) * t194 + (r_i_i_C(1) + pkin(8) + pkin(7) + qJ(2)) * qJD(1) + t186 - t212;
t195 = sin(qJ(1));
t216 = t194 * t195;
t211 = t188 * t216;
t196 = cos(qJ(1));
t214 = qJD(1) * t196;
t227 = t187 * t214 + t211;
t222 = pkin(3) * t189;
t218 = t188 * t194;
t217 = t194 * t187;
t215 = qJD(1) * t195;
t213 = qJD(5) * t188;
t210 = t196 * t218;
t207 = t187 * t215;
t209 = pkin(4) * t207 + r_i_i_C(2) * t210 + t196 * t213;
t205 = t220 * t187;
t203 = t227 * r_i_i_C(2) + t195 * t213 + t204 * t214;
t202 = -r_i_i_C(2) * t187 - t204;
t200 = -t223 * t217 + t220 * t218 + t186;
t199 = (-pkin(4) * t188 - t205) * t194;
t190 = cos(t193);
t198 = -t190 * t219 + t199;
t197 = qJD(2) + (-t223 * t188 - pkin(3) * t190 - cos(pkin(10)) * pkin(2) - pkin(1) - t205) * qJD(1);
t1 = [-t195 * t228 + t197 * t196, t214, t198 * t196 + (t202 + t222) * t215 + t209, t196 * t199 + t202 * t215 + t209, -t207 + t210, 0; t197 * t195 + t196 * t228, t215 (-pkin(4) * t187 - t222) * t214 + t198 * t195 + t203, -pkin(4) * t211 + (-pkin(4) * t214 - t220 * t216) * t187 + t203, t227, 0; 0, 0, t200 - t212, t200, t217, 0;];
JaD_transl  = t1;
