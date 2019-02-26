% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR14
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:45
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRPRPR14_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR14_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR14_jacobig_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:45:26
% EndTime: 2019-02-26 21:45:26
% DurationCPUTime: 0.03s
% Computational Cost: add. (8->8), mult. (24->18), div. (0->0), fcn. (40->8), ass. (0->15)
t101 = sin(pkin(6));
t105 = sin(qJ(1));
t114 = t105 * t101;
t104 = sin(qJ(2));
t113 = t105 * t104;
t107 = cos(qJ(2));
t112 = t105 * t107;
t108 = cos(qJ(1));
t111 = t108 * t101;
t110 = t108 * t104;
t109 = t108 * t107;
t106 = cos(qJ(4));
t103 = sin(qJ(4));
t102 = cos(pkin(6));
t1 = [0, t114, 0, -t102 * t113 + t109, 0, t106 * t114 + (t102 * t112 + t110) * t103; 0, -t111, 0, t102 * t110 + t112, 0, -t106 * t111 + (-t102 * t109 + t113) * t103; 1, t102, 0, t101 * t104, 0, -t101 * t107 * t103 + t102 * t106;];
Jg_rot  = t1;
