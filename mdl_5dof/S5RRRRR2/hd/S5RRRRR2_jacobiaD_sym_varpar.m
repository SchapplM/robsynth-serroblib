% Zeitableitung der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRRRR2
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S5RRRRR2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% 
% Output:
% JaD [6x5]
%   Zeitableitung der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-29 15:26
% Revision: 932832b1be1be80f59b7f1a581a1a8f328bdb39d (2019-03-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD = S5RRRRR2_jacobiaD_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)



%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR2_jacobiaD_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR2_jacobiaD_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRRR2_jacobiaD_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRR2_jacobiaD_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5RRRRR2_jacobiaD_sym_varpar: pkin has to be [2x1] (double)');
%% Function calls
if link_index == 0
	JaD=S5RRRRR2_jacobiaD_0_sym_varpar(qJ, qJD, r_i_i_C, pkin);
elseif link_index == 1
	JaD=S5RRRRR2_jacobiaD_1_sym_varpar(qJ, qJD, r_i_i_C, pkin);
elseif link_index == 2
	JaD=S5RRRRR2_jacobiaD_2_sym_varpar(qJ, qJD, r_i_i_C, pkin);
elseif link_index == 3
	JaD=S5RRRRR2_jacobiaD_3_sym_varpar(qJ, qJD, r_i_i_C, pkin);
elseif link_index == 4
	JaD=S5RRRRR2_jacobiaD_4_sym_varpar(qJ, qJD, r_i_i_C, pkin);
elseif link_index == 5
	JaD=S5RRRRR2_jacobiaD_5_sym_varpar(qJ, qJD, r_i_i_C, pkin);
else
	JaD=NaN(6,5);
end