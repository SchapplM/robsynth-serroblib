% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPPRPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:29
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPRPR7_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR7_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR7_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRPR7_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR7_jacobiaD_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:29:10
% EndTime: 2019-02-26 20:29:10
% DurationCPUTime: 0.15s
% Computational Cost: add. (121->30), mult. (196->45), div. (0->0), fcn. (147->7), ass. (0->21)
t171 = pkin(9) + qJ(4);
t168 = sin(t171);
t169 = cos(t171);
t172 = sin(pkin(10));
t174 = cos(pkin(10));
t181 = r_i_i_C(1) * t174 - r_i_i_C(2) * t172 + pkin(4);
t188 = r_i_i_C(3) + qJ(5);
t192 = -(t188 * t168 + t181 * t169) * qJD(4) + t169 * qJD(5);
t191 = t181 * qJD(4) - qJD(5);
t189 = t181 * t168 - t188 * t169 + sin(pkin(9)) * pkin(3) + qJ(2);
t176 = sin(qJ(1));
t187 = qJD(1) * t176;
t177 = cos(qJ(1));
t170 = qJD(1) * t177;
t186 = qJD(4) * t176;
t185 = qJD(4) * t177;
t182 = qJD(1) * t188;
t180 = -t172 * r_i_i_C(1) - t174 * r_i_i_C(2) - pkin(1) - pkin(7) - qJ(3);
t179 = qJD(1) * t181;
t178 = qJD(2) - t192;
t1 = [-t176 * qJD(3) + t178 * t177 + (-t189 * t176 + t180 * t177) * qJD(1), t170, -t187 (t177 * t179 + t188 * t186) * t169 + (-t191 * t176 + t177 * t182) * t168, t168 * t186 - t169 * t170, 0; t177 * qJD(3) + t178 * t176 + (t180 * t176 + t189 * t177) * qJD(1), t187, t170 (t176 * t179 - t188 * t185) * t169 + (t176 * t182 + t191 * t177) * t168, -t168 * t185 - t169 * t187, 0; 0, 0, 0, t192, qJD(4) * t169, 0;];
JaD_transl  = t1;
