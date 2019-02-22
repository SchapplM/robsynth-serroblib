% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPRP4
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 10:27
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRPRP4_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP4_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP4_jacobiR_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 10:27:36
% EndTime: 2019-02-22 10:27:36
% DurationCPUTime: 0.06s
% Computational Cost: add. (36->11), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->16)
t67 = sin(qJ(5));
t68 = sin(qJ(3));
t74 = t68 * t67;
t69 = cos(qJ(5));
t73 = t68 * t69;
t70 = cos(qJ(3));
t72 = t70 * t67;
t71 = t70 * t69;
t66 = qJ(1) + pkin(9);
t65 = cos(t66);
t64 = sin(t66);
t63 = -t64 * t74 + t65 * t69;
t62 = t64 * t73 + t65 * t67;
t61 = t64 * t69 + t65 * t74;
t60 = -t64 * t67 + t65 * t73;
t1 = [t63, 0, t65 * t72, 0, t60, 0; t61, 0, t64 * t72, 0, t62, 0; 0, 0, t74, 0, -t71, 0; -t62, 0, t65 * t71, 0, -t61, 0; t60, 0, t64 * t71, 0, t63, 0; 0, 0, t73, 0, t72, 0; -t64 * t70, 0, -t65 * t68, 0, 0, 0; t65 * t70, 0, -t64 * t68, 0, 0, 0; 0, 0, t70, 0, 0, 0;];
JR_rot  = t1;
