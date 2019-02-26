% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:33
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRPPRR9_jacobig_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR9_jacobig_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR9_jacobig_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:33:12
% EndTime: 2019-02-26 21:33:12
% DurationCPUTime: 0.02s
% Computational Cost: add. (3->3), mult. (9->8), div. (0->0), fcn. (17->6), ass. (0->8)
t72 = cos(pkin(6));
t75 = cos(qJ(2));
t77 = t72 * t75;
t76 = cos(qJ(1));
t74 = sin(qJ(1));
t73 = sin(qJ(2));
t71 = sin(pkin(6));
t1 = [0, t74 * t71, 0, 0, -t76 * t73 - t74 * t77, 0; 0, -t76 * t71, 0, 0, -t74 * t73 + t76 * t77, 0; 1, t72, 0, 0, t71 * t75, 0;];
Jg_rot  = t1;
