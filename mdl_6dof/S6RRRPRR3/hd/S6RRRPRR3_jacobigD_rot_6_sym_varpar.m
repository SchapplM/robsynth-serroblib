% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR3
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:17
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRPRR3_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR3_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR3_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR3_jacobigD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:17:08
% EndTime: 2019-02-26 22:17:08
% DurationCPUTime: 0.07s
% Computational Cost: add. (37->8), mult. (50->12), div. (0->0), fcn. (50->6), ass. (0->13)
t181 = qJD(5) - qJD(2) - qJD(3);
t175 = sin(qJ(1));
t168 = qJD(1) * t175;
t177 = cos(qJ(1));
t169 = qJD(1) * t177;
t173 = qJ(2) + qJ(3);
t170 = sin(t173);
t171 = cos(t173);
t174 = sin(qJ(5));
t176 = cos(qJ(5));
t180 = t170 * t176 - t171 * t174;
t178 = t181 * (t170 * t174 + t171 * t176);
t1 = [0, t169, t169, 0, -t169, t180 * t168 + t178 * t177; 0, t168, t168, 0, -t168, -t180 * t169 + t178 * t175; 0, 0, 0, 0, 0, t181 * t180;];
JgD_rot  = t1;
