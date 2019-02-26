% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRRR4
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:20
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRRRR4_jacobigD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR4_jacobigD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR4_jacobigD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR4_jacobigD_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:20:27
% EndTime: 2019-02-26 20:20:27
% DurationCPUTime: 0.06s
% Computational Cost: add. (44->17), mult. (161->43), div. (0->0), fcn. (169->10), ass. (0->27)
t209 = sin(pkin(7));
t210 = sin(pkin(6));
t228 = t210 * t209;
t212 = cos(pkin(7));
t216 = cos(qJ(3));
t227 = t212 * t216;
t213 = cos(pkin(6));
t215 = sin(qJ(2));
t226 = t213 * t215;
t217 = cos(qJ(2));
t225 = t213 * t217;
t214 = sin(qJ(3));
t224 = t214 * t217;
t223 = t215 * t216;
t222 = qJD(2) * t214;
t208 = sin(pkin(13));
t211 = cos(pkin(13));
t221 = -t208 * t215 + t211 * t225;
t220 = t208 * t217 + t211 * t226;
t219 = -t208 * t225 - t211 * t215;
t218 = t208 * t226 - t211 * t217;
t207 = t218 * qJD(2);
t206 = t220 * qJD(2);
t205 = t213 * t209 * qJD(3) * t214 + ((t212 * t224 + t223) * qJD(3) + (t212 * t223 + t224) * qJD(2)) * t210;
t204 = -t207 * t227 + t219 * t222 + (-t218 * t216 + (t208 * t228 + t219 * t212) * t214) * qJD(3);
t203 = t206 * t227 + t221 * t222 + (t220 * t216 + (-t211 * t228 + t221 * t212) * t214) * qJD(3);
t1 = [0, 0, -t207 * t209, t204, t204, 0; 0, 0, t206 * t209, t203, t203, 0; 0, 0, qJD(2) * t215 * t228, t205, t205, 0;];
JgD_rot  = t1;
