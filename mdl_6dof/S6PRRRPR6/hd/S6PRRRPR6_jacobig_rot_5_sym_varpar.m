% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRPR6
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:13
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRRRPR6_jacobig_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR6_jacobig_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR6_jacobig_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:13:25
% EndTime: 2019-02-26 20:13:25
% DurationCPUTime: 0.03s
% Computational Cost: add. (9->9), mult. (24->20), div. (0->0), fcn. (40->8), ass. (0->13)
t111 = sin(pkin(11));
t112 = sin(pkin(6));
t122 = t111 * t112;
t113 = cos(pkin(11));
t121 = t113 * t112;
t114 = cos(pkin(6));
t116 = sin(qJ(2));
t120 = t114 * t116;
t118 = cos(qJ(2));
t119 = t114 * t118;
t117 = cos(qJ(3));
t115 = sin(qJ(3));
t1 = [0, t122, t111 * t119 + t113 * t116 (-t111 * t120 + t113 * t118) * t115 - t117 * t122, 0, 0; 0, -t121, t111 * t116 - t113 * t119 (t111 * t118 + t113 * t120) * t115 + t117 * t121, 0, 0; 0, t114, -t112 * t118, t112 * t116 * t115 - t114 * t117, 0, 0;];
Jg_rot  = t1;
