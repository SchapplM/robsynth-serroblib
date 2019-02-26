% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPPRRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:30
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPRRP2_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP2_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP2_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRP2_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP2_jacobiaD_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:30:42
% EndTime: 2019-02-26 20:30:42
% DurationCPUTime: 0.17s
% Computational Cost: add. (246->45), mult. (280->75), div. (0->0), fcn. (219->9), ass. (0->34)
t217 = pkin(10) + qJ(4);
t213 = sin(t217);
t215 = cos(t217);
t240 = pkin(8) + r_i_i_C(3);
t232 = t240 * t215;
t244 = (-pkin(4) * t213 + t232) * qJD(4);
t221 = cos(qJ(5));
t229 = qJD(1) * t215 - qJD(5);
t242 = t221 * t229;
t233 = qJD(5) * t215;
t230 = -qJD(1) + t233;
t220 = sin(qJ(5));
t236 = qJD(4) * t220;
t241 = -t213 * t236 + t230 * t221;
t218 = qJ(1) + pkin(9);
t214 = sin(t218);
t238 = qJD(1) * t214;
t216 = cos(t218);
t237 = qJD(1) * t216;
t235 = qJD(4) * t221;
t234 = qJD(5) * t213;
t228 = r_i_i_C(1) * t220 + r_i_i_C(2) * t221;
t227 = r_i_i_C(1) * t221 - r_i_i_C(2) * t220 + pkin(4);
t226 = t229 * t220;
t225 = -pkin(4) * t215 - t240 * t213 - cos(pkin(10)) * pkin(3) - pkin(2);
t224 = qJD(4) * t227;
t223 = t213 * t235 + t230 * t220;
t222 = -t240 * qJD(4) + t228 * qJD(5);
t219 = -pkin(7) - qJ(3);
t211 = t223 * t214 - t216 * t242;
t210 = t241 * t214 + t216 * t226;
t209 = t214 * t242 + t223 * t216;
t208 = t214 * t226 - t241 * t216;
t1 = [t211 * r_i_i_C(1) + t210 * r_i_i_C(2) + t216 * qJD(3) - t214 * t244 + (-cos(qJ(1)) * pkin(1) + t214 * t219 + t225 * t216) * qJD(1), 0, t237 (-t216 * t224 - t240 * t238) * t215 + (t222 * t216 + t227 * t238) * t213, t208 * r_i_i_C(1) + t209 * r_i_i_C(2), 0; -t209 * r_i_i_C(1) + t208 * r_i_i_C(2) + t214 * qJD(3) + t216 * t244 + (-sin(qJ(1)) * pkin(1) - t216 * t219 + t225 * t214) * qJD(1), 0, t238 (-t214 * t224 + t240 * t237) * t215 + (t222 * t214 - t227 * t237) * t213, -t210 * r_i_i_C(1) + t211 * r_i_i_C(2), 0; 0, 0, 0, -t228 * t233 + (-t227 * t213 + t232) * qJD(4) (-t215 * t235 + t220 * t234) * r_i_i_C(2) + (-t215 * t236 - t221 * t234) * r_i_i_C(1), 0;];
JaD_transl  = t1;
