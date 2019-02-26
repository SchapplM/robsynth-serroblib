% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRPR1
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:46
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRPRPR1_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR1_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_jacobigD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:46:28
% EndTime: 2019-02-26 19:46:28
% DurationCPUTime: 0.06s
% Computational Cost: add. (36->15), mult. (97->37), div. (0->0), fcn. (104->10), ass. (0->20)
t185 = sin(pkin(11));
t188 = cos(pkin(11));
t191 = sin(qJ(2));
t192 = cos(qJ(2));
t194 = t191 * t185 - t192 * t188;
t197 = t194 * qJD(2);
t195 = t185 * t192 + t188 * t191;
t180 = t195 * qJD(2);
t184 = qJ(4) + pkin(12);
t182 = sin(t184);
t187 = sin(pkin(6));
t196 = t187 * t182;
t190 = cos(pkin(6));
t189 = cos(pkin(10));
t186 = sin(pkin(10));
t183 = cos(t184);
t178 = t195 * t190;
t177 = t190 * t197;
t176 = t190 * t180;
t1 = [0, 0, 0, -t186 * t176 - t189 * t197, 0 (t186 * t177 - t189 * t180) * t182 + ((-t186 * t178 - t189 * t194) * t183 + t186 * t196) * qJD(4); 0, 0, 0, t189 * t176 - t186 * t197, 0 (-t189 * t177 - t186 * t180) * t182 + ((t189 * t178 - t186 * t194) * t183 - t189 * t196) * qJD(4); 0, 0, 0, t187 * t180, 0, t190 * qJD(4) * t182 + (t195 * qJD(4) * t183 - t182 * t197) * t187;];
JgD_rot  = t1;
