% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPP7
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:59
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRPP7_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP7_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP7_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPP7_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP7_jacobiaD_transl_6_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:59:52
% EndTime: 2019-02-26 20:59:52
% DurationCPUTime: 0.32s
% Computational Cost: add. (215->56), mult. (654->85), div. (0->0), fcn. (546->6), ass. (0->40)
t209 = sin(qJ(3));
t212 = cos(qJ(3));
t208 = sin(qJ(4));
t211 = cos(qJ(4));
t230 = r_i_i_C(1) + pkin(5) + pkin(4);
t234 = qJD(5) * t208;
t244 = r_i_i_C(2) + qJ(5);
t216 = -t234 + (t230 * t208 - t244 * t211) * qJD(4);
t217 = t244 * t208 + t230 * t211 + pkin(3);
t229 = pkin(8) - r_i_i_C(3) - qJ(6);
t246 = t229 * t212;
t253 = (t217 * t209 - t246) * qJD(3) + qJD(6) * t209 + t216 * t212;
t222 = t229 * t209;
t250 = t217 * t212 + t222;
t249 = (-pkin(3) * t209 - qJ(2) + t246) * qJD(1) + t211 * qJD(5);
t210 = sin(qJ(1));
t243 = t210 * t208;
t242 = t210 * t211;
t213 = cos(qJ(1));
t241 = t213 * t208;
t240 = qJD(1) * t210;
t239 = qJD(1) * t213;
t238 = qJD(3) * t209;
t237 = qJD(3) * t212;
t236 = qJD(3) * t213;
t235 = qJD(4) * t212;
t231 = t212 * qJD(6);
t228 = t213 * t209 * t211;
t227 = t212 * t236;
t225 = qJD(4) * t209 + qJD(1);
t224 = qJD(1) * t209 + qJD(4);
t220 = t224 * t213;
t219 = t209 * t242 + t241;
t215 = qJD(1) * t250;
t214 = t209 * t234 + t231 + qJD(2) + (-pkin(1) - pkin(7)) * qJD(1) + (pkin(3) * t212 + t222) * qJD(3);
t203 = t211 * t220 + (-t225 * t208 + t211 * t237) * t210;
t202 = t225 * t242 + (t210 * t237 + t220) * t208;
t201 = -t211 * t227 + (t209 * t241 + t242) * qJD(4) + t219 * qJD(1);
t200 = -qJD(4) * t228 - t208 * t227 - t211 * t239 + t224 * t243;
t1 = [-t244 * t200 - t230 * t201 + t249 * t210 + t214 * t213, t239, -t253 * t210 + t213 * t215, t219 * qJD(5) - t230 * t202 + t244 * t203, t202, -t210 * t238 + t212 * t239; t244 * t202 + t230 * t203 + t214 * t210 - t249 * t213, t240, t210 * t215 + t253 * t213 -(t228 - t243) * qJD(5) + t244 * t201 - t230 * t200, t200, t209 * t236 + t212 * t240; 0, 0, -t250 * qJD(3) + t216 * t209 - t231 (t230 * t238 - t244 * t235) * t208 + (-t244 * t238 + (-t230 * qJD(4) + qJD(5)) * t212) * t211, -t208 * t238 + t211 * t235, -t237;];
JaD_transl  = t1;
